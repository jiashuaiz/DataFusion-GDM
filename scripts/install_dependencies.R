# Install required R packages for DataFusion-GDM

cran_core <- c(
  "ggplot2",
  "vegan",
  "mice",
  "VIM"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Installing %s...", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

invisible(lapply(cran_core, install_if_missing))

if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  message("Installing optional heatmap dependencies via Bioconductor (ComplexHeatmap, circlize) and RColorBrewer...")
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install(c("ComplexHeatmap", "circlize"), ask = FALSE, update = FALSE)
  install_if_missing("RColorBrewer")
}

message("All core dependencies checked/installed. Optional heatmap deps installed if requested.")
