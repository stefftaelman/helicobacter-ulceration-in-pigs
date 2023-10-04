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
metadata$fg <- ifelse(metadata$FG_PO_QPCR > 0, 1, 0)
dim(metadata)
dim(raw_counts)

relative_counts <- prop.table(as.matrix(raw_counts), 2)

# create phyloseq object
ps <- phyloseq::phyloseq(
  otu_table=phyloseq::otu_table(raw_counts, taxa_are_rows=TRUE),
  sample_data=phyloseq::sample_data(metadata),
  tax_table=phyloseq::tax_table(taxonomy)
)

# agglomerate taxa
ps_rel <- phyloseq::transform_sample_counts(ps, function(x) x/sum(x))
genus_glom <- phyloseq::tax_glom(ps_rel, taxrank = 'genus')
# create dataframe from phyloseq object
genus_abundances <- data.table::data.table(psmelt(genus_glom))
# convert genus to a character vector from a factor because R
genus_abundances$genus <- as.character(genus_abundances$genus)
# group dataframe by genus, calculate median rel. abundance
genus_abundances[, mean := mean(Abundance, na.rm = TRUE), by = "genus"]
# Change name of lower abundant genera
genus_abundances[(mean <= 0.014), genus := "Other"]
unique(genus_abundances$genus)
# filter and sort on abundance
genus_abundances <- genus_abundances %>%
  filter(Abundance > 0)
sorted_genera <- genus_abundances %>%
  dplyr::count(genus, wt = Abundance) %>%
  arrange(desc(n)) %>%
  pull(genus)
# put "Other" last in sorted genera
sorted_genera <- c(sorted_genera[sorted_genera != "Other"], "Other")
sorted_samples <- genus_abundances %>%
  filter(genus == sorted_genera[1]) %>%
  arrange (desc(Abundance)) %>%
  pull(Sample) %>%
  unique(.)
genus_names <- factor(
  c("Lactobacillus", "Turicibacter", "Streptococcus", "Clostridium sensu stricto 1", "Unnamed_Clostridiaceae", "Unnamed_Clostridia", "Terrisporobacter", "Unnamed_Peptostreptococcaceae", "Actinobacillus", "Unnamed_Pasteurellaceae", "Escherichia-Shigella", "Moraxella", "Fusobacterium", "Unnamed_Bacteroidales", "Unnamed_organism", "Other"),
  levels=c("Lactobacillus", "Turicibacter", "Streptococcus", "Clostridium sensu stricto 1", "Unnamed_Clostridiaceae", "Unnamed_Clostridia", "Terrisporobacter", "Unnamed_Peptostreptococcaceae", "Actinobacillus", "Unnamed_Pasteurellaceae", "Escherichia-Shigella", "Moraxella", "Fusobacterium", "Unnamed_Bacteroidales", "Unnamed_organism", "Other")
)
genus_colors <- c("#00A98F", "#80D297", "#98dfb5", "#F7E733", "#FFC000", "#FF9F1C", "#1565A9", "#8FBFE0", "#C00000", "#F0544F", "#7E52A0", "#FFC6D9", "#D4F5F5", "#4C3B4D", "#0D0D0D", "#E1E1E1")
genus_abundances <- genus_abundances %>%
  mutate(Sample = factor(Sample, levels = sorted_samples)) %>%
  mutate(genus = factor(genus, levels = rev(genus_names)))
tax_barplot <- genus_abundances %>%
  ggplot2::ggplot(ggplot2::aes(x=Sample, y=Abundance)) +
  ggplot2::geom_bar(
    ggplot2::aes(fill=genus),
    stat = 'identity',
    position = 'fill'
    ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 60, size = 5, vjust = 1, hjust = 1
        ),
    legend.position = 'bottom',
    legend.title = ggplot2::element_text(size=12),
    legend.key.size = ggplot2::unit(0.5, 'cm'),
    legend.text = ggplot2::element_text(size = 10)
    ) + 
  ggplot2::labs(x="", y="Relative abundance") +
  ggplot2::scale_fill_manual(
    "Genus",
    breaks=genus_names,
    values=genus_colors,
    labels=c(
      "Lactobacillus", "Turicibacter", "Streptococcus",
      "Clostridium sensu stricto 1", "Unnamed Clostridiaceae", "Unnamed Clostridia",
      "Terrisporobacter", "Unnamed Peptostreptococcaceae",
      "Actinobacillus", "Unnamed Pasteurellaceae",
      "Escherichia-Shigella", "Moraxella",
      "Fusobacterium", "Unnamed Bacteroidales",
      "Unnamed organism", "Other"
    )
  ) + 
  ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0))
tax_barplot

svglite::svglite(
  filename = "./figures/taxonomy_barplot_genus.svg",
  width=10, height=8
  )
tax_barplot
dev.off()