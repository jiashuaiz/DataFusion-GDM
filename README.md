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
# Optional: install dependencies first
source("scripts/install_dependencies.R")

# Install from local checkout
install.packages(".", repos = NULL, type = "source")

# Or using devtools
# devtools::install_local(".")
```

## Usage

```r
library(DataFusionGDM)

# Simulate a GDM and save to CSV
export_simulated_gdm("GDM_simulated.csv", scenario = "default", n_pops = 40)

# Run a scenario and get plots
res <- run_genetic_scenario("island", n_pops = 40)
res$plots$heatmap()
res$plots$mds()
```

## Examples

Reproducible analysis pipelines are provided under `inst/examples`:
- `01_mds_procrustes_sensitivity.R` — MDS + Procrustes sensitivity study
- `02_besmi_prepare_and_batch.R` — BESMI dataset preparation and batch imputation

To run after installation:

```r
example_file <- system.file("examples", "01_mds_procrustes_sensitivity.R", package = "DataFusionGDM")
source(example_file)
```

## Contents

- Simulation and visualization APIs in `R/simulate_gdm.R`
- MDS & Procrustes APIs in `R/mds_procrustes.R`
- BESMI preparation and imputation APIs in `R/besmi*.R`
- Example data in `inst/extdata/GDM_simulated.csv`

## License

GPL-3.0
