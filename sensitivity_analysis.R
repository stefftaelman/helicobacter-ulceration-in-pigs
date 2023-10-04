library(dplyr)

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
alpha_diversity$stable <- as.factor(alpha_diversity$stable)
row.names(alpha_diversity) <- alpha_diversity$Row.names
alpha_diversity$Row.names <- NULL
head(alpha_diversity)

# association between hs in the fundus and shannon diversity
m_fhs_sha <- glm(
  Shannon ~ fundus_hs + stable,
  data=alpha_diversity,
  family="gaussian"
  )
summary(m_fhs_sha)

sensitivity <- sensemakr::sensemakr(
    model = m_fhs_sha,
    treatment = "fundus_hs",
    benchmark_covariates = "stable21",
    kd = c(8, 8.5, 9)
)

png(
  filename = "./figures/fundusHS_shannondiv_sensitivity_analysis.png",
  width=800, height=800, res=150
  )
plot(sensitivity)

dev.off()

# association between hs in the pylorus and shannon diversity
m_phs_sha <- glm(
  Shannon ~ pylorus_hs + stable,
  data=alpha_diversity,
  family="gaussian"
  )
summary(m_phs_sha)

sensitivity <- sensemakr::sensemakr(
    model = m_phs_sha,
    treatment = "pylorus_hs",
    benchmark_covariates = "stable21",
    kd = c(35, 40, 45)
)

png(
  filename = "./figures/pylorusHS_shannondiv_sensitivity_analysis.png",
  width=800, height=800, res=150
  )
plot(sensitivity)
dev.off()
