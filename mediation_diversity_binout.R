library(dplyr)
library(sl3)
library(medoutcon)
set.seed(42)

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
taxonomy_path <- "./data/asv_taxonomy.csv"

# read in data
metadata <- read.csv(metadata_path, row.names = 1)
raw_counts <- data.frame(t(data.table::fread(raw_counts_path)[, -1]))
seqs <- row.names(raw_counts)
row.names(raw_counts) <- paste0("ASV", seq(dim(raw_counts)[1]))
colnames(raw_counts) <- row.names(metadata)
taxonomy <- as.matrix(read.csv(taxonomy_path, row.names = 1))[seqs, ]
row.names(taxonomy) <- row.names(raw_counts)

internal_controls <- row.names(metadata[metadata$sample_type == "control", ])
print(sum(raw_counts[, internal_controls]))
raw_counts <- raw_counts[, !names(raw_counts) %in% internal_controls]
metadata <- metadata[!rownames(metadata) %in% internal_controls, !names(metadata) %in% c("fastq_1", "fastq_2", "sample_type")]
metadata$stable <- as.factor(metadata$stable)
dim(metadata)
dim(raw_counts)

## calculate alpha diversity
ps_object <- phyloseq::phyloseq(
  otu_table=phyloseq::otu_table(raw_counts, taxa_are_rows=TRUE),
  sample_data=phyloseq::sample_data(metadata),
  tax_table=phyloseq::tax_table(taxonomy)
)
alpha_diversity <- phyloseq::estimate_richness(ps_object, split=TRUE, "Shannon")
alpha_diversity <- merge(x=alpha_diversity, y=metadata, by="row.names")
alpha_diversity <- alpha_diversity %>%
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
row.names(alpha_diversity) <- alpha_diversity$Row.names
alpha_diversity$Row.names <- NULL
head(alpha_diversity)

# Mediation analysis PGZ
## for categorical outcome variable
lrn_ridge_cat <- sl3::Lrnr_glmnet$new(alpha = 0, family = "multinomial")
lrn_lasso_cat <- sl3::Lrnr_glmnet$new(alpha = 1, family = "multinomial")
stack_cat <- sl3::Stack$new(lrn_ridge_cat, lrn_lasso_cat)
sl_cat <- sl3::Lrnr_sl$new(learners = stack_cat)
## for binary output variable
lrn_glm_bin <- sl3::Lrnr_glm$new(family = "binomial")
lrn_ridge_bin <- sl3::Lrnr_glmnet$new(alpha = 0, family = "binomial")
lrn_lasso_bin <- sl3::Lrnr_glmnet$new(alpha = 1, family = "binomial")
lrn_ranger <- sl3::Lrnr_ranger$new()
stack_bin <- sl3::Stack$new(lrn_glm_bin, lrn_ridge_bin, lrn_lasso_bin, lrn_ranger)
sl_bin <- sl3::Lrnr_sl$new(learners = stack_bin, metalearner = Lrnr_nnls$new())
## for continuous output variable
lrn_ridge_contin <- sl3::Lrnr_glmnet$new(alpha=0, family="gaussian")
lrn_lasso_contin <- sl3::Lrnr_glmnet$new(alpha=1, family="gaussian")
lrn_fglm_contin <- sl3::Lrnr_glm_fast$new(family=gaussian())
lrn_ranger <- sl3::Lrnr_ranger$new()
stack_contin <- sl3::Stack$new(lrn_ridge_contin, lrn_lasso_contin, lrn_fglm_contin, lrn_ranger)
sl_contin <- sl3::Lrnr_sl$new(learners = stack_contin, metalearner = Lrnr_nnls$new())

w <- cbind(
  #alpha_diversity$stable,
  alpha_diversity$stable11,alpha_diversity$stable12,alpha_diversity$stable13,alpha_diversity$stable14,alpha_diversity$stable21,alpha_diversity$stable23,
  alpha_diversity$pellet_feed
  )
a <- alpha_diversity$pylorus_hs
y <- alpha_diversity$ulcer_binary
m <- alpha_diversity$Shannon

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

# Mediation analysis FGZ
w <- cbind(
  #alpha_diversity$stable,
  alpha_diversity$stable11,alpha_diversity$stable12,alpha_diversity$stable13,alpha_diversity$stable14,alpha_diversity$stable21,alpha_diversity$stable23,
  alpha_diversity$pellet_feed
  )
a <- alpha_diversity$fundus_hs
y <- alpha_diversity$ulcer_binary
m <- alpha_diversity$Shannon

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
