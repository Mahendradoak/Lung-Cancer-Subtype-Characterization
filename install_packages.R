# Installs everything the pipeline needs. Tested with R 4.3.3 / Bioconductor 3.18:
# limma 3.58.1, ConsensusClusterPlus 1.66.0, survival 3.5-8, cluster 2.1.6,
# data.table, matrixStats, jsonlite, testthat 3.2.1, withr.
# MOFA2 (training factors) also needs Python via basilisk; clusterProfiler and
# org.Hs.eg.db are optional (GO enrichment is skipped without them).
cran <- c("data.table", "matrixStats", "jsonlite", "survival", "cluster", "testthat", "withr")
install.packages(setdiff(cran, rownames(installed.packages())), repos = "https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
bioc <- c("limma", "ConsensusClusterPlus", "MOFA2", "clusterProfiler", "org.Hs.eg.db")
BiocManager::install(setdiff(bioc, rownames(installed.packages())), update = FALSE, ask = FALSE)
