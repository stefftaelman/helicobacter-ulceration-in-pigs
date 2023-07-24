library(dplyr)
library(sl3)
library(medoutcon)

COLOR_SCHEME <- c(
  "#00A98F",
  "#FFC000",
  "#1565A9",
  "#C00000",
  "#0D0D0D"
)

# set paths
metadata_path <- "./data/metadata.csv"
raw_counts_path <- "./data/asv_counts.csv"
taxonomy_path <- "./data/taxonomy_clean.csv"

# read in data
metadata <- read.csv(metadata_path, row.names = 1)
raw_counts <- data.frame(t(data.table::fread(raw_counts_path)[, -1]))
seqs <- row.names(raw_counts)
row.names(raw_counts) <- paste0("ASV", seq(dim(raw_counts)[1]))
colnames(raw_counts) <- row.names(metadata)
taxonomy <- as.matrix(read.csv(taxonomy_path, row.names = 1))[row.names(raw_counts), ]

internal_controls <- row.names(metadata[metadata$sample_type == "control", ])
print(sum(raw_counts[, internal_controls]))
raw_counts <- raw_counts[, !names(raw_counts) %in% internal_controls]
metadata <- metadata[!rownames(metadata) %in% internal_controls, !names(metadata) %in% c("fastq_1", "fastq_2", "sample_type")]
metadata$stable <- as.factor(metadata$stable)
dim(metadata)
dim(raw_counts)

# read in DAA taxa
daa_taxa_pgz <- read.csv("./deliverables/diff_abundance_score_adjusted_tmm_pgz.csv", row.names = 1)
key_taxa_pgz <- daa_taxa_pgz[daa_taxa_pgz$passes_permFDP == TRUE, ]$names
daa_taxa_fgz <- read.csv("./deliverables/diff_abundance_score_adjusted_tmm_fgz.csv", row.names = 1)
key_taxa_fgz <- daa_taxa_fgz[daa_taxa_fgz$passes_permFDP == TRUE, ]$names
key_taxa_names <- intersect(key_taxa_pgz, key_taxa_fgz)
key_taxa_IDs <- row.names(taxonomy[taxonomy[, ncol(taxonomy)] %in% key_taxa_names, ])
enriched_taxa_names <- daa_taxa_pgz[(daa_taxa_pgz$logFC >= 2 & daa_taxa_pgz$names %in% key_taxa_names), ]$names
enriched_taxa_IDs <- row.names(taxonomy[taxonomy[, ncol(taxonomy)] %in% enriched_taxa_names, ])
diminished_taxa_names <- daa_taxa_pgz[(daa_taxa_pgz$logFC <= -2 & daa_taxa_pgz$names %in% key_taxa_names), ]$names
diminished_taxa_IDs <- row.names(taxonomy[taxonomy[, ncol(taxonomy)] %in% diminished_taxa_names, ])

key_taxa_counts <- t(raw_counts[key_taxa_IDs, ])

## extract key species diversities
metadata <- merge(x=metadata, y=key_taxa_counts, by="row.names")
metadata <- metadata %>%
    dplyr::mutate(
        pylorus_hs=ifelse(HS_PGZ_QPCR > 0, 1, 0),
        fundus_hs=ifelse(HS_FGZ_QPCR > 0, 1, 0),
        pellet_feed=ifelse(group == "Korrel", 1, 0),
        stable11=ifelse(stable == 11, 1, 0),
        stable12=ifelse(stable == 12, 1, 0),
        stable13=ifelse(stable == 13, 1, 0),
        stable14=ifelse(stable == 14, 1, 0),
        stable21=ifelse(stable == 21, 1, 0),
        stable23=ifelse(stable == 23, 1, 0),
        ulcer_binary=ifelse(ulcer_score > 2, 1, 0)
        )
row.names(metadata) <- metadata$Row.names
metadata$Row.names <- NULL
head(metadata)

# Enriched species mediation
## Mediation analysis PGZ
### for categorical outcome variable
lrn_ridge_cat <- sl3::Lrnr_glmnet$new(alpha = 0, family = "multinomial")
lrn_lasso_cat <- sl3::Lrnr_glmnet$new(alpha = 1, family = "multinomial")
stack_cat <- sl3::Stack$new(lrn_ridge_cat, lrn_lasso_cat)
sl_cat <- sl3::Lrnr_sl$new(learners = stack_cat)
### for binary output variable
lrn_glm_bin <- sl3::Lrnr_glm$new(family = "binomial")
lrn_ridge_bin <- sl3::Lrnr_glmnet$new(alpha = 0, family = "binomial")
lrn_lasso_bin <- sl3::Lrnr_glmnet$new(alpha = 1, family = "binomial")
lrn_ranger <- sl3::Lrnr_ranger$new()
stack_bin <- sl3::Stack$new(lrn_glm_bin, lrn_ridge_bin, lrn_lasso_bin, lrn_ranger)
sl_bin <- sl3::Lrnr_sl$new(learners = stack_bin, metalearner = Lrnr_nnls$new())
### for continuous output variable
lrn_ridge_contin <- sl3::Lrnr_glmnet$new(alpha=0, family="gaussian")
lrn_lasso_contin <- sl3::Lrnr_glmnet$new(alpha=1, family="gaussian")
lrn_fglm_contin <- sl3::Lrnr_glm_fast$new(family=gaussian())
lrn_ranger <- sl3::Lrnr_ranger$new()
stack_contin <- sl3::Stack$new(lrn_ridge_contin, lrn_lasso_contin, lrn_fglm_contin, lrn_ranger)
sl_contin <- sl3::Lrnr_sl$new(learners = stack_contin, metalearner = Lrnr_nnls$new())

w <- cbind(
  #metadata$stable,
  metadata$stable11,metadata$stable12,metadata$stable13,metadata$stable14,metadata$stable21,metadata$stable23,
  metadata$pellet_feed
  )
a <- metadata$pylorus_hs
y <- metadata$ulcer_binary
m <- metadata[, enriched_taxa_IDs]

tmle_dne_pgz <- medoutcon(
  W = w, A = a, Z = NULL, M = m, Y = y, 
  effect = "direct", 
  estimator = "tmle", 
  g_learners = sl_bin,
  h_learners = sl_bin,
  b_learners = sl_bin
  )
tmle_dne_pgz

tmle_ine_pgz <- medoutcon(
  W = w, A = a, Z = NULL, M = m, Y = y, 
  effect = "indirect", 
  estimator = "tmle", 
  g_learners = sl_bin,
  h_learners = sl_bin,
  b_learners = sl_bin
  )
tmle_ine_pgz

## Mediation analysis FGZ
w <- cbind(
  #metadata$stable,
  metadata$stable11,metadata$stable12,metadata$stable13,metadata$stable14,metadata$stable21,metadata$stable23,
  metadata$pellet_feed
  )
a <- metadata$fundus_hs
y <- metadata$ulcer_binary
m <- metadata[, enriched_taxa_IDs]

tmle_dne_fgz <- medoutcon(
  W = w, A = a, Z = NULL, M = m, Y = y, 
  effect = "direct", 
  estimator = "tmle", 
  g_learners = sl_bin,
  h_learners = sl_bin,
  b_learners = sl_bin
  )
tmle_dne_fgz

tmle_ine_fgz <- medoutcon(
  W = w, A = a, Z = NULL, M = m, Y = y, 
  effect = "indirect", 
  estimator = "tmle", 
  g_learners = sl_bin,
  h_learners = sl_bin,
  b_learners = sl_bin
  )
tmle_ine_fgz

# Diminished species mediation
## Mediation analysis PGZ
w <- cbind(
  #metadata$stable,
  metadata$stable11,metadata$stable12,metadata$stable13,metadata$stable14,metadata$stable21,metadata$stable23,
  metadata$pellet_feed
  )
a <- metadata$pylorus_hs
y <- metadata$ulcer_binary
m <- metadata[, dimished_taxa_IDs]

tmle_dne_pgz <- medoutcon(
  W = w, A = a, Z = NULL, M = m, Y = y, 
  effect = "direct", 
  estimator = "tmle", 
  g_learners = sl_bin,
  h_learners = sl_bin,
  b_learners = sl_bin
  )
tmle_dne_pgz

tmle_ine_pgz <- medoutcon(
  W = w, A = a, Z = NULL, M = m, Y = y, 
  effect = "indirect", 
  estimator = "tmle", 
  g_learners = sl_bin,
  h_learners = sl_bin,
  b_learners = sl_bin
  )
tmle_ine_pgz

## Mediation analysis FGZ
w <- cbind(
  #metadata$stable,
  metadata$stable11,metadata$stable12,metadata$stable13,metadata$stable14,metadata$stable21,metadata$stable23,
  metadata$pellet_feed
  )
a <- metadata$fundus_hs
y <- metadata$ulcer_binary
m <- metadata[, dimished_taxa_IDs]

tmle_dne_fgz <- medoutcon(
  W = w, A = a, Z = NULL, M = m, Y = y, 
  effect = "direct", 
  estimator = "tmle", 
  g_learners = sl_bin,
  h_learners = sl_bin,
  b_learners = sl_bin
  )
tmle_dne_fgz

tmle_ine_fgz <- medoutcon(
  W = w, A = a, Z = NULL, M = m, Y = y, 
  effect = "indirect", 
  estimator = "tmle", 
  g_learners = sl_bin,
  h_learners = sl_bin,
  b_learners = sl_bin
  )
tmle_ine_fgz
