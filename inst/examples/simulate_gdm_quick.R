library(DataFusionGDM)

# Quick simulation and visualization
res <- run_genetic_scenario("island", n_pops = 40, seed = 123)
print(res$plots$heatmap())
print(res$plots$mds())

