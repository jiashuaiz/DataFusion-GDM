# DataFusion-GDM

- Purpose: Machine Learning Solutions for Integrating Partially Overlapped Genetic Datasets
- Version: V1.0.0
- Author: Zhu, Jiashuai (The University of Melbourne; Agriculture Victoria Research)

## Genetic Distance Simulation (`00-Simulate_structured_distM`)

Creates a realistic genetic distance matrix (GDM) with assumed population structure:
- Generates population metadata with hierarchical group structures
- Places populations in multidimensional space based on genetic and geographic factors
- Supports creation of admixed populations and bottleneck effects
- Calculates distances between populations with realistic biological properties
- Outputs a complete GDM to serve as ground truth for imputation evaluation

### Execution
- run `Simulate_structured_distM.r`

## MDS and Procrustes Analysis (`01-MDS_Procrustes`)

A pipeline for aligning and comparing distance matrices through dimensional reduction:
- Performs Multidimensional Scaling (MDS) on distance matrices
- Applies Procrustes transformation to align MDS spaces
- Includes sensitivity analysis to evaluate numbers of shared populations needed
- Visualizes the sensitivity analysis
- Converts transformed MDSs back to distance matrices

### Execution

1. Ensure a complete genetic distance matrix (`GDM_simulated.csv`), if using controlled population structures, `input` 
2. Run `03_procrustes_sensitivity.r`
	- `03_procrustes_sensitivity.r` will call `01_prepare_data_xxx.r` and `02_perform_mds+procrustes.r`
3. Run `04_visualization.r` to present results

## Bootstrapping Evaluation for Structural Missingness Imputation (BESMI) (`02-BESMI`)

A framework for evaluating GDM imputation with structured missingness:
- Creates test datasets with controlled patterns of missing values
- Implements multiple imputation functions with convergence tracking
- Processes datasets in batch mode with checkpointing

## Execution Flow

1. Ensure a genetic distance matrix (`data/GDM_simulated.csv`)
2. Run `01-BESMI_Prepare_datasets.r` to create test datasets
3. Run `03-Batch_Processing.r` to perform imputation across all datasets; this will call `02-Imputation.r` in batch
