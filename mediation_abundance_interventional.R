library(dplyr)
library(sl3)
library(doFuture)
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
dim(metadata)
dim(raw_counts)

# read in DAA taxa
# daa_taxa_pgz <- read.csv("./deliverables/diff_abundance_score_adjusted_tmm_pgz.csv", row.names = 1)
# key_taxa_pgz <- daa_taxa_pgz[daa_taxa_pgz$passes_permFDP == TRUE, ]$names
# daa_taxa_fgz <- read.csv("./deliverables/diff_abundance_score_adjusted_tmm_fgz.csv", row.names = 1)
# key_taxa_fgz <- daa_taxa_fgz[daa_taxa_fgz$passes_permFDP == TRUE, ]$names
# key_taxa_names <- intersect(key_taxa_pgz, key_taxa_fgz)
# key_taxa_IDs <- row.names(taxonomy[taxonomy[, ncol(taxonomy)] %in% key_taxa_names, ])
# enriched_taxa_names <- daa_taxa_pgz[(daa_taxa_pgz$logFC >= 2 & daa_taxa_pgz$names %in% key_taxa_names), ]$names
# enriched_taxa_IDs <- names(which(taxonomy[, ncol(taxonomy)] == enriched_taxa_names[1]))
# diminished_taxa_names <- daa_taxa_pgz[(daa_taxa_pgz$logFC <= -2 & daa_taxa_pgz$names %in% key_taxa_names), ]$names
# diminished_taxa_IDs <- row.names(taxonomy[taxonomy[, ncol(taxonomy)] %in% diminished_taxa_names, ])


ps <- phyloseq::phyloseq(
  phyloseq::otu_table(raw_counts, taxa_are_rows=TRUE),
  phyloseq::sample_data(metadata),
  phyloseq::tax_table(taxonomy)
  )

## Mediation analysis PGZ
extract_vars <- function(ps, gland_zone="PGZ", taxonomic_level="genus"){
  stopifnot(taxonomic_level %in% phyloseq::rank_names(ps))
  stopifnot(gland_zone %in% c("PGZ", "FGZ"))
  # Extract the exposure
  if (gland_zone == "PGZ"){
    A <- phyloseq::sample_data(ps)$pylorus_hs
  } else if (gland_zone == "FGZ"){
    A <- phyloseq::sample_data(ps)$fundus_hs
  }
  # Extract the outcome
  Y <- phyloseq::sample_data(ps)$ulcer_binary
  # Extract the mediators
  tax_level_ps <- phyloseq::tax_glom(ps, taxonomic_level)
  M <- t(phyloseq::otu_table(tax_level_ps))
  colnames(M) <- phyloseq::tax_table(tax_level_ps)[, taxonomic_level]
  # prevalence filtering (cutoff 90% => change to 5%)
  M <- as.matrix(M[, apply(M, 2, function(col) mean(col != 0)) >= 0.05])
  # extract the covariates
  covars <- c(
    "stable11", "stable12", "stable13", "stable14", "stable21", "stable23"
    )
  W <- data.frame(phyloseq::sample_data(ps)) %>%
    select(all_of(covars))
  if (impute){
    W <- W %>%
      mutate(across(everything(), impute_lowest))
  }
  # take complete cases
  complete_cases <- complete.cases(A, Y, M, W)
  if (sum(complete_cases) != length(A)){
    message("Removed ", length(A) - sum(complete_cases), " observations due to missing values.")
  }
  A <- A[complete_cases]
  Y <- Y[complete_cases]
  M <- M[complete_cases, ]
  W <- W[complete_cases, ]
  # Convert covariates to numeric
  W <- data.frame(lapply(W, function(x) as.numeric(as.character(x)))) %>%
    as.matrix()
  message(
    paste0(
      "Extracted exposure A with ", length(A), " observations.\n",
      "Extracted outcome Y with ", length(Y), " observations.\n",
      "Extracted mediators M with ", dim(M)[1], " observations and ", dim(M)[2], " features.\n",
      "Extracted covariates W with ", dim(W)[1], " observations and ", dim(W)[2], " features."
      )
    )
  vars_df <- data.frame(A, Y, M, W)
  return(
    list(
      df=vars_df,
      mediator_names=names(vars_df)[3:(ncol(vars_df)-ncol(W))],
      covariate_names=names(vars_df)[(ncol(vars_df)-ncol(W)+1):ncol(vars_df)]
    )
  )
}

#sl_library <- c("SL.glm","SL.step","SL.glm.interaction","SL.randomForest","SL.mean")
sl_library <- c("SL.glm","SL.mean")
vars_pgz <- extract_vars(ps, gland_zone="PGZ", taxonomic_level="genus")

mediation_results_pgz <- data.frame()
for (i in seq_along(vars_pgz$mediator_names)){
  ans <- HDmediation::mediation(
    data = vars_pgz$df,
    A = "A", W = vars_pgz$covariate_names,
    Z = vars_pgz$mediator_names[-i],
    M = vars_pgz$mediator_names[i],
    Y = "Y", S = NULL,
    family = "gaussian", folds = 1, partial_tmle = TRUE, bounds = NULL,
    learners_g = sl_library, learners_e = sl_library, learners_c = sl_library,
    learners_b = sl_library, learners_hz = sl_library, learners_u = sl_library,
    learners_ubar = sl_library, learners_v = sl_library, learners_vbar = sl_library
  )
  ans$mediator <- vars_pgz$mediator_names[i]
  mediation_results_pgz <- rbind(mediation_results_pgz, ans)
}
write.csv(mediation_results_pgz, "./deliverables/taxon_mediation_results_pgz.csv")

## Mediation analysis FGZ
vars_fgz <- extract_vars(ps, gland_zone="FGZ", taxonomic_level="genus")

mediation_results_fgz <- data.frame()
for (i in seq_along(vars_fgz$mediator_names)){
  ans <- HDmediation::mediation(
    data = vars_fgz$df,
    A = "A", W = vars_fgz$covariate_names,
    Z = vars_fgz$mediator_names[-i],
    M = vars_fgz$mediator_names[i],
    Y = "Y", S = NULL,
    family = "gaussian", folds = 1, partial_tmle = TRUE, bounds = NULL,
    learners_g = sl_library, learners_e = sl_library, learners_c = sl_library,
    learners_b = sl_library, learners_hz = sl_library, learners_u = sl_library,
    learners_ubar = sl_library, learners_v = sl_library, learners_vbar = sl_library
  )
  ans$mediator <- vars_fgz$mediator_names[i]
  mediation_results_fgz <- rbind(mediation_results_fgz, ans)
}
write.csv(mediation_results_fgz, "./deliverables/taxon_mediation_results_fgz.csv")