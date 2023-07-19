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
metadata$ulcer_score <- factor(metadata$ulcer_score, levels=c("2", "3", "4", "5"))
metadata$stable <- as.factor(metadata$stable)
dim(metadata)
dim(raw_counts)


### binary (qualitative) regressions
hs_fgz_bin <- ifelse(metadata$HS_FGZ_QPCR > 0, 1, 0)
hs_pgz_bin <- ifelse(metadata$HS_PGZ_QPCR > 0, 1, 0)
fg_po_bin <- ifelse(metadata$FG_PO_QPCR > 0, 1, 0)
pellet_feed <- ifelse(metadata$group == "Korrel", 1, 0)

# single effect regressions
## association between diet and ulceration score
m_diet_ulcer <- MASS::polr(factor(metadata$ulcer_score) ~ pellet_feed, method="logistic")
summary(m_diet_ulcer)

## association between diet and abundance of hs in the fundus
m_diet_hs_fgz <- glm(metadata$HS_FGZ_QPCR ~ pellet_feed, family="gaussian")
summary(m_diet_hs_fgz)

## association between diet and abundance of hs in the pylorus
m_diet_hs_pgz <- glm(metadata$HS_PGZ_QPCR ~ pellet_feed, family="gaussian")
summary(m_diet_hs_pgz)

## association between hs and fundus fg
m_hs_fg <- glm(metadata$FG_PO_QPCR ~ metadata$HS_FGZ_QPCR + factor(metadata$stable), family="gaussian")
summary(m_hs_fg)

## association between hs and pylorus fg
m_hs_fg <- glm(metadata$FG_PO_QPCR ~ metadata$HS_PGZ_QPCR + factor(metadata$stable), family="gaussian")
summary(m_hs_fg)

# association with microbiome diversity
## calculate Shannon diversity
ps_object <- phyloseq::phyloseq(
  otu_table=phyloseq::otu_table(raw_counts, taxa_are_rows=TRUE),
  sample_data=phyloseq::sample_data(metadata),
  tax_table=phyloseq::tax_table(taxonomy)
)
alpha_diversity <- phyloseq::estimate_richness(ps_object, split=TRUE, "Shannon")
alpha_diversity <- merge(x=alpha_diversity, y=metadata, by="row.names")
alpha_diversity <- alpha_diversity %>%
    dplyr::mutate(
        pylorus_hs=ifelse(HS_PGZ_QPCR > 0, "Present", "Absent"),
        fundus_hs=ifelse(HS_FGZ_QPCR > 0, "Present", "Absent"),
        text=paste0(
            "Sample: ", Row.names, " (stable: ", stable, ")\n",
            "Shannon diversity index: ", round(Shannon, digits=3), "\n",
            "H. Suis qPCR count: ", round(HS_FGZ_QPCR, digits=3), " (fundus); ", round(HS_PGZ_QPCR, digits=3), " (pylorus)\n",
            "F. Gastrosuis qPCR count: ", FG_PO_QPCR, "\n"
        )
    )
row.names(alpha_diversity) <- alpha_diversity$Row.names
alpha_diversity$Row.names <- NULL

boxplot_pgz <- ggplot2::ggplot(
    data=alpha_diversity,
    ggplot2::aes(x=pylorus_hs, y=Shannon)
    ) +
    ggplot2::geom_boxplot(alpha=0.2) +
    ggplot2::geom_jitter(
        alpha=0.9,
        ggplot2::aes(shape=stable, color=pylorus_hs, size=HS_PGZ_QPCR, text=text),
    ) +
    ggplot2::theme_bw() +
    ggplot2::coord_flip() +
    ggplot2::scale_color_manual(values=COLOR_SCHEME[1:2]) +
    ggplot2::labs(x="", y="Shannon diversity index")

png("./figures/diversity_HS_PGZ_presence_boxplot.png", width=800, height=600, res=100)
boxplot_pgz
dev.off()

svglite::svglite("./figures/diversity_HS_PGZ_presence_boxplot.svg", width=8, height=6)
boxplot_pgz
dev.off()

boxplot_fgz <- ggplot2::ggplot(
    data=alpha_diversity,
    ggplot2::aes(x=fundus_hs, y=Shannon)
    ) +
    ggplot2::geom_boxplot(alpha=0.2) +
    ggplot2::geom_jitter(
        alpha=0.9,
        ggplot2::aes(shape=stable, color=fundus_hs, size=HS_FGZ_QPCR, text=text),
    ) +
    ggplot2::theme_bw() +
    ggplot2::coord_flip() +
    ggplot2::scale_color_manual(values=COLOR_SCHEME[1:2]) +
    ggplot2::labs(x="", y="Shannon diversity index")

png("./figures/diversity_HS_FGZ_presence_boxplot.png", width=800, height=600, res=100)
boxplot_fgz
dev.off()

svglite::svglite("./figures/diversity_HS_FGZ_presence_boxplot.svg", width=8, height=6)
boxplot_fgz
dev.off()

## association between hs and shannon diversity
m_div_hs <- glm(alpha_diversity$Shannon ~ alpha_diversity$HS_FGZ_QPCR + factor(alpha_diversity$stable), family="gaussian")
summary(m_div_hs)

m_div_hs <- glm(alpha_diversity$Shannon ~ alpha_diversity$HS_PGZ_QPCR + factor(alpha_diversity$stable), family="gaussian")
summary(m_div_hs)

m_div_hs_bin <- glm(alpha_diversity$Shannon ~ ifelse(alpha_diversity$HS_FGZ_QPCR > 0, 1, 0) + factor(alpha_diversity$stable), family="gaussian")
summary(m_div_hs_bin)

m_div_hs_bin <- glm(metadalpha_diversityata$Shannon ~ ifelse(alpha_diversity$HS_PGZ_QPCR > 0, 1, 0) + factor(alpha_diversity$stable), family="gaussian")
summary(m_div_hs_bin)
