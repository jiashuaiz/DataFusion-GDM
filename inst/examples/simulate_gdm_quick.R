library(DataFusionGDM)

# Quick simulation and visualization
set.seed(123)
res <- run_genetic_scenario("island", n_pops = 40)
print(res$plots$heatmap())
print(res$plots$mds())

