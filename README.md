# DataFusionGDM  
[![DOI (software)](https://img.shields.io/badge/DOI-10.26188/28602953-blue)](https://doi.org/10.26188/28602953)
[![DOI (paper)](https://img.shields.io/badge/DOI-10.3389%2Ffpls.2025.1543956-blue)](https://doi.org/10.3389/fpls.2025.1543956)

Machine Learning Solutions for Integrating Partially Overlapped Genetic Datasets.

## How to cite

If you use DataFusion-GDM, please cite:

- Paper: Zhu J., Malmberg M.M., Shinozuka M., Retegan R.M., Cogan N.O., Jacobs J.L., Giri K., Smith K.F. (2025). Machine learning solutions for integrating partially overlapping genetic datasets and modelling host–endophyte effects in ryegrass (Lolium) dry matter yield estimation. Frontiers in Plant Science. https://doi.org/10.3389/fpls.2025.1543956
- Software: Zhu, J. (2025). DataFusion-GDM. The University of Melbourne. Software. https://doi.org/10.26188/28602953

Author ORCID: https://orcid.org/0000-0002-9916-9732

## Install

In R:

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("jiashuaiz/DataFusion-GDM")
```

## Usage

```r
library(DataFusionGDM)

# Simulate a GDM in memory and visualize
res <- run_genetic_scenario("island", n_pops = 40)
res$plots$heatmap()
res$plots$mds()

# Optionally export to CSV if needed (defaults to tempdir)
tmp <- export_simulated_gdm(scenario = "default", n_pops = 40, verbose = FALSE)
# unlink(tmp)  # clean up when finished

# Simulate and visualize
source(system.file("examples/simulate_gdm_quick.R", package = "DataFusionGDM"), echo = TRUE)

# MDS + Procrustes
source(system.file("examples/mds_procrustes_demo.R", package = "DataFusionGDM"), echo = TRUE)

# BESMI batch (small demo)
source(system.file("examples/besmi_batch_quick.R", package = "DataFusionGDM"), echo = TRUE)
```

## Vignettes

See the package vignettes for end-to-end guides:
- Getting started
- MDS + Procrustes sensitivity
- BESMI batch imputation

Open vignettes in R:

```r
browseVignettes("DataFusionGDM")
vignette("getting-started", package = "DataFusionGDM")
```

## Contents

- Simulation and visualization APIs in `R/simulate_gdm.R`
- MDS & Procrustes APIs in `R/mds_procrustes.R`
- BESMI preparation and imputation APIs in `R/besmi*.R`
 - Vignettes under `vignettes/` (no bundled data; examples use in-memory/temp files)

## License

GPL-3.0
