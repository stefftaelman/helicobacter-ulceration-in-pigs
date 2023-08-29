library(dplyr)

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
metadata$ulcer_score <- factor(metadata$ulcer_score, levels=c("2", "3", "4", "5"))
metadata$stable <- as.factor(metadata$stable)
metadata$pylorus_hs <- ifelse(metadata$HS_PGZ_QPCR > 0, 1, 0)
metadata$fundus_hs <- ifelse(metadata$HS_FGZ_QPCR > 0, 1, 0)
dim(metadata)
dim(raw_counts)

phylogenetic_tree <- TreeSummarizedExperiment::toTree(
  data = taxonomy
  )

# Differential abundance analysis
row.names(taxonomy) <- taxonomy[, ncol(taxonomy)]
row.names(raw_counts) <- row.names(taxonomy)

## TSE object
tse_object <- TreeSummarizedExperiment::TreeSummarizedExperiment(
  assays = list(count = as.matrix(raw_counts)),
  rowData = taxonomy,
  colData = metadata,
  rowTree = phylogenetic_tree,
  rowNodeLab = rownames(raw_counts)
)

## Run DA PGZ
res_score_pgz <- treeclimbR::runDA(
  TSE = tse_object,
  feature_on_row = TRUE,
  assay = 1,
  #design_terms = "ulcer_score"
  design_terms = c("stable", "pylorus_hs", "ulcer_score")
  )
res_table_score_pgz <- treeclimbR::nodeResult(
  object = res_score_pgz,
  n = Inf
  )
candidates_score_pgz <- treeclimbR::getCand(
  tree = TreeSummarizedExperiment::rowTree(tse_object),
  score_data = res_table_score_pgz,
  node_column = "node",
  p_column = "PValue",
  sign_column = "logFC"
  )
da_results_score_pgz <- treeclimbR::evalCand(
  tree = TreeSummarizedExperiment::rowTree(tse_object),
  type = "single",
  levels = candidates_score_pgz$candidate_list,
  limit_rej = 0.1,
  score_data = res_table_score_pgz,
  node_column = "node",
  p_column = "PValue",
  sign_column = "logFC",
  use_pseudo_leaf = FALSE
  )

da_results_score_pgz$output$names <- TreeSummarizedExperiment::convertNode(
  phylogenetic_tree,
  da_results_score_pgz$output$node
  )


### permutation-based FDR estimation
fdr_thresh_score_pgz <- permFDP::permFDP.adjust.threshold(
  pVals = da_results_score_pgz$output$PValue,
  threshold = 0.05,
  myDesign = ifelse(metadata[colnames(raw_counts), "ulcer_score"] == 2, 1, 2),
  intOnly = raw_counts[da_results_score_pgz$output$names, ],
  nPerms = 999
)
da_results_score_pgz$output$passes_permFDP <- da_results_score_pgz$output$PValue <= fdr_thresh_score_pgz
write.csv(da_results_score_pgz$output, "./deliverables/diff_abundance_score_adjusted_tmm_pgz.csv")
DT::datatable(data.frame(da_results_score_pgz$output))

### volcanoplot

#' Create a volcano plot from the outputs of a differential abundance analysis
#'
#' @param results A list containing the DA results as outputted in the
#' `treeclimbR::evalCand` function.
#' @param alpha a numeric alpha level on which to cut off p-values
#' @param effect_cutoff a numeric log fold change on which to cut off log FCs
#' @param color_scheme a list of hex valued colors indicating a scheme.
#' @param title title for the plot to export
#'
#' @examples volcano_plot(results, color_scheme)
#' @export
volcano_plot <- function(results,
                         alpha = 0.05,
                         effect_cutoff = 2,
                         color_scheme, 
                         title = "Differential abundance analysis: Volcano plot") {
  # extract elements
  scores <- results$output

  # color on significance groups and add tooltips
  volcano_plot_df <- scores %>%
    mutate(
      significance_group = dplyr::if_else(
        abs(logFC) < effect_cutoff & PValue > alpha,
        "Non-significant",
        dplyr::if_else(
          abs(logFC) < effect_cutoff,
          "Statistically significant (p-values)",
          dplyr::if_else(
            PValue > alpha,
            "Biologically significant (FC)",
            "Biologically & statistically significant"
          )
        )
      ),
      text = paste0(
        "Taxon: ", names, "\n",
        "p-value: ", PValue, "\n",
        "log FC: ", round(logFC, 4), "\n"
      )
    )

  p <- ggplot2::ggplot(
    volcano_plot_df,
    ggplot2::aes(
      x = logFC, y = -log10(PValue),
      color = significance_group,
      text = text
    )
  ) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(
      values = c(
        "Non-significant" = "#9EA0A5",
        "Biologically significant (FC)" = color_scheme[1],
        "Statistically significant (p-values)" = color_scheme[3],
        "Biologically & statistically significant" = color_scheme[2]
      )
    ) +
    ggplot2::labs(
      title = title,
      x = "Fold change (Log2)",
      y = "p-value (-Log10)"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(vjust = 0.4, hjust = 1),
      axis.text.y = ggplot2::element_text(vjust = 0.4, hjust = 1),
      legend.title = ggplot2::element_blank(),
      legend.position = "right",
      text = ggplot2::element_text(size = 10)
    ) +
    ggplot2::geom_vline(
      xintercept = c(-effect_cutoff, effect_cutoff),
      linetype = "dashed"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(alpha),
      linetype = "dashed"
    )

  return(plotly::ggplotly(p, tooltip = "text"))
}

p_score_pgz <- volcano_plot(
  da_results_score_pgz,
  alpha = fdr_thresh_score_pgz, 
  color_scheme = COLOR_SCHEME, 
  title = "Differential taxa associated to ulcer score"
  )

png("./figures/DAA_ulceration_adjusted_pgz_volcanoplot.png", width=800, height=600)
p_score_pgz
dev.off()

svglite::svglite("./figures/DAA_ulceration_adjusted_pgz_volcanoplot.svg", width=8, height=6)
p_score_pgz
dev.off()

htmlwidgets::saveWidget(
  p_score_pgz,
  file = "./figures/DAA_ulceration_adjusted_pgz_volcanoplot.html",
  selfcontained=F
  )

## Run DA FGZ
res_score_fgz <- treeclimbR::runDA(
  TSE = tse_object,
  feature_on_row = TRUE,
  assay = 1,
  #design_terms = "ulcer_score"
  design_terms = c("stable", "fundus_hs", "ulcer_score")
  )
res_table_score_fgz <- treeclimbR::nodeResult(
  object = res_score_fgz,
  n = Inf
  )
candidates_score_fgz <- treeclimbR::getCand(
  tree = TreeSummarizedExperiment::rowTree(tse_object),
  score_data = res_table_score_fgz,
  node_column = "node",
  p_column = "PValue",
  sign_column = "logFC"
  )
da_results_score_fgz <- treeclimbR::evalCand(
  tree = TreeSummarizedExperiment::rowTree(tse_object),
  type = "single",
  levels = candidates_score_fgz$candidate_list,
  limit_rej = 0.1,
  score_data = res_table_score_fgz,
  node_column = "node",
  p_column = "PValue",
  sign_column = "logFC",
  use_pseudo_leaf = FALSE
  )

da_results_score_fgz$output$names <- TreeSummarizedExperiment::convertNode(
  phylogenetic_tree,
  da_results_score_fgz$output$node
  )


### permutation-based FDR estimation
fdr_thresh_score_fgz <- permFDP::permFDP.adjust.threshold(
  pVals = da_results_score_fgz$output$PValue,
  threshold = 0.05,
  myDesign = ifelse(metadata[colnames(raw_counts), "ulcer_score"] == 2, 1, 2),
  intOnly = raw_counts[da_results_score_fgz$output$names, ],
  nPerms = 999
)
da_results_score_fgz$output$passes_permFDP <- da_results_score_fgz$output$PValue <= fdr_thresh_score_fgz
write.csv(da_results_score_fgz$output, "./deliverables/diff_abundance_score_adjusted_tmm_fgz.csv")
DT::datatable(data.frame(da_results_score_fgz$output))

### volcanoplot
p_score_fgz <- volcano_plot(
  da_results_score_fgz,
  alpha = fdr_thresh_score_fgz, 
  color_scheme = COLOR_SCHEME, 
  title = "Differential taxa associated to ulcer score"
  )

png("./figures/DAA_ulceration_adjusted_fgz_volcanoplot.png", width=800, height=600)
p_score_fgz
dev.off()

svglite::svglite("./figures/DAA_ulceration_adjusted_fgz_volcanoplot.svg", width=8, height=6)
p_score_fgz
dev.off()

htmlwidgets::saveWidget(
  p_score_fgz,
  file = "./figures/DAA_ulceration_adjusted_fgz_volcanoplot.html",
  selfcontained=F
  )
