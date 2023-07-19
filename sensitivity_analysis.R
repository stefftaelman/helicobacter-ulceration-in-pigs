metadata_path <- "./data/metadata.csv"

metadata <- read.csv(metadata_path, row.names = 1)
metadata <- metadata[complete.cases(metadata), ]
head(metadata)


### binary (qualitative) regressions
hs_fgz_bin <- ifelse(metadata$HS_FGZ_QPCR > 0, 1, 0)
hs_pgz_bin <- ifelse(metadata$HS_PGZ_QPCR > 0, 1, 0)
fg_po_bin <- ifelse(metadata$FG_PO_QPCR > 0, 1, 0)
pellet_feed <- ifelse(metadata$group == "Korrel", 1, 0)

# association between hs and fg
m_hs_fg <- glm(FG_PO_QPCR ~ HS_FGZ_QPCR + stable, data=metadata ,family="gaussian")
summary(m_hs_fg)

sensitivity <- sensemakr::sensemakr(
    model = m_hs_fg,
    treatment = "HS_FGZ_QPCR",
    benchmark_covariates = "stable",
    kd = c(5, 10)
)

png(
  filename = "./figures/HS_FG_sensitivity_analysis.png",
  width=800, height=800, res=150
  )
plot(sensitivity)
dev.off()
