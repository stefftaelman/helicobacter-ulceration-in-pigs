library(data.table)
library(dplyr)
library(sl3)
library(doFuture)
library(ranger)
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

ps <- phyloseq::phyloseq(
  phyloseq::otu_table(raw_counts, taxa_are_rows=TRUE),
  phyloseq::sample_data(metadata),
  phyloseq::tax_table(taxonomy)
  )

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

find_significant_effects <- function(mediation_results, top_n=NULL){
  # look for mediators that have confidence intervals (for either indirect or for direct) that do not include zero
  significant_effects <- mediation_results %>%
    filter(
      (0 %between% list(ci_indirect_low, ci_indirect_high)) == FALSE #|
      #  (0 %between% list(ci_direct_low, ci_direct_high)) == FALSE
      )
  # sort by the biggest difference between indirect and direct effect
  significant_effects <- significant_effects %>%
    mutate(diff = abs(indirect - direct)) %>%
    arrange(desc(diff))
  if (!is.null(top_n)){
    significant_effects <- significant_effects %>%
      slice(1:top_n)
  }
  return(significant_effects)
}

plot_significant_effects <- function(mediation_results, color_scheme=COLOR_SCHEME, top_n=NULL, export_name=NULL){
  significant_effects <- find_significant_effects(mediation_results, top_n=top_n)

  # Melt the dataframe to long format
  long_effects <- reshape2::melt(
    significant_effects,
    id.vars = "mediator",
    measure.vars = c("indirect", "direct"),
    variable.name = "effect_type",
    value.name = "effect"
  )
  indirect_ci_data <- data.frame(
    mediator = significant_effects$mediator,
    effect_type = "indirect",
    lower_ci = significant_effects$ci_indirect_low,
    upper_ci = significant_effects$ci_indirect_high
  )
  direct_ci_data <- data.frame(
    mediator = significant_effects$mediator,
    effect_type = "direct",
    lower_ci = significant_effects$ci_direct_low,
    upper_ci = significant_effects$ci_direct_high
  )
  ci_data <- rbind(indirect_ci_data, direct_ci_data)
  # merge the long effects dataframe with the confidence intervals dataframe
  long_effects <- merge(long_effects, ci_data, by = c("mediator", "effect_type"))

  # Create the plot
  p <- ggplot2::ggplot(
    long_effects,
    ggplot2::aes(x = effect, y = reorder(mediator, abs(effect)), color = effect_type)
    ) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = lower_ci, xmax = upper_ci),
      width = 0.2, position = ggplot2::position_dodge(width = 0.5)
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::labs(
      x = "Effect estimate",
      y = "Mediator",
      color = "Effect Type"
      ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_manual(
      name = "",values = color_scheme[c(2, 1)]
      )
  if (!is.na(export_name)) {
    ggplot2::ggsave(export_name, p)
  }
  return(p)
}

## Mediation analysis PGZ

# check if the taxon_mediation_results_pgz.csv file exists
if (!file.exists("./deliverables/taxon_mediation_results_pgz.csv")){
  sl_library <- c("SL.glm","SL.mean", "SL.ranger")
  vars_pgz <- extract_vars(ps, gland_zone="PGZ", taxonomic_level="genus")

  mediation_results_pgz <- data.frame()
  for (i in seq_along(vars_pgz$mediator_names)){
    cat(paste0(round(i/length(vars_pgz$mediator_names) * 100), '% completed'))
    ans <- HDmediation::mediation(
      data = vars_pgz$df,
      A = "A", W = vars_pgz$covariate_names,
      Z = vars_pgz$mediator_names[-i],
      M = vars_pgz$mediator_names[i],
      Y = "Y", S = NULL,
      family = "binomial", folds = 2, partial_tmle = TRUE, bounds = NULL,
      learners_g = sl_library, learners_e = sl_library, learners_c = sl_library,
      learners_b = sl_library, learners_hz = sl_library, learners_u = sl_library,
      learners_ubar = sl_library, learners_v = sl_library, learners_vbar = sl_library
    )
    ans$mediator <- vars_pgz$mediator_names[i]
    mediation_results_pgz <- rbind(mediation_results_pgz, ans)
    if (i == length(vars_pgz$mediator_names)) cat(': Done')
    else cat('\014')
  }
  write.csv(mediation_results_pgz, "./deliverables/taxon_mediation_results_pgz.csv")
} else {
  mediation_results_pgz <- read.csv("./deliverables/taxon_mediation_results_pgz.csv", row.names = 1, header=TRUE)
}
p_pgz <- plot_significant_effects(mediation_results_pgz, top_n=20, export_name="./deliverables/mediation_effects_pgz.svg")
plotly::ggplotly(p_pgz)

## Mediation analysis FGZ
# check if the taxon_mediation_results_fgz.csv file exists
if (!file.exists("./deliverables/taxon_mediation_results_fgz.csv")){
  vars_fgz <- extract_vars(ps, gland_zone="FGZ", taxonomic_level="genus")

  mediation_results_fgz <- data.frame()
  for (i in seq_along(vars_fgz$mediator_names)){
    cat(paste0(round(i/length(vars_fgz$mediator_names) * 100), '% completed'))
    ans <- HDmediation::mediation(
      data = vars_fgz$df,
      A = "A", W = vars_fgz$covariate_names,
      Z = vars_fgz$mediator_names[-i],
      M = vars_fgz$mediator_names[i],
      Y = "Y", S = NULL,
      family = "binomial", folds = 1, partial_tmle = TRUE, bounds = NULL,
      learners_g = sl_library, learners_e = sl_library, learners_c = sl_library,
      learners_b = sl_library, learners_hz = sl_library, learners_u = sl_library,
      learners_ubar = sl_library, learners_v = sl_library, learners_vbar = sl_library
    )
    ans$mediator <- vars_fgz$mediator_names[i]
    mediation_results_fgz <- rbind(mediation_results_fgz, ans)
    if (i == length(vars_fgz$mediator_names)) cat(': Done')
    else cat('\014')
  }
  write.csv(mediation_results_fgz, "./deliverables/taxon_mediation_results_fgz.csv")
} else {
  mediation_results_fgz <- read.csv("./deliverables/taxon_mediation_results_fgz.csv", row.names = 1, header=TRUE)
}

p_fgz <- plot_significant_effects(mediation_results_fgz, top_n=15, export_name="./deliverables/mediation_effects_fgz.svg")
plotly::ggplotly(p_fgz)
