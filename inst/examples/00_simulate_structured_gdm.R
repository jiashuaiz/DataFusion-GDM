# Example: Structured genetic distance simulation (concise)

cat("[DataFusionGDM] Running genetics-only scenario...\n")
library(DataFusionGDM)

set.seed(233)
res <- run_genetic_scenario("genetics", n_pops = 30)

cat("[DataFusionGDM] Showing heatmap and MDS...\n")
res$plots$heatmap()
res$plots$mds()

cat("[DataFusionGDM] Done.\n")


