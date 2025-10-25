# DataFusionGDM

Machine Learning Solutions for Integrating Partially Overlapped Genetic Datasets.

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
