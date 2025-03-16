# 04_visualization.r
# Purpose: Visualize sensitivity analysis results for different shared populations
# V1.0.0
# Author: Zhu, Jiashuai (The University of Melbourne; Agriculture Victoria Research)

rm(list = ls())

# Load required libraries
library(ggplot2)
library(scales)

data_source <- "input"  # Options: "random" or "input"

# Configure parameters based on data source
if (data_source == "random") {
  results_df <- readRDS("data/sensitivity_summary_random.rds")
} else if (data_source == "input") {
  results_df <- readRDS("data/sensitivity_summary_input.rds")
} else {
  stop("Invalid data_source. Choose 'random' or 'input'")
}

# Normalize k to get proportion of shared populations
results_df$prop <- results_df$k / max(results_df$k)

# Define plot limits
xlim_min <- min(results_df$prop)

# Plot the Prior and Posterior
p1 <- ggplot(results_df, aes(x = prop)) +
  geom_line(aes(y = the_prior.mean, color = "Prior"), linewidth = 1) +
  geom_point(aes(y = the_prior.mean, color = "Prior"), size = 2) +
  geom_ribbon(aes(ymin = the_prior.mean - the_prior.sd, 
                  ymax = the_prior.mean + the_prior.sd, 
                  fill = "Prior"), alpha = 0.2) +
  
  geom_line(aes(y = the_posterior.mean, color = "Posterior"), linewidth = 1) +
  geom_point(aes(y = the_posterior.mean, color = "Posterior"), size = 2) +
  geom_ribbon(aes(ymin = the_posterior.mean - the_posterior.sd, 
                  ymax = the_posterior.mean + the_posterior.sd, 
                  fill = "Posterior"), alpha = 0.2) +
  
  scale_x_continuous(labels = percent, limits = c(xlim_min, 1)) +
  labs(title = "the Prior vs Posterior",
       x = "Proportion of Shared Populations",
       y = "Differences") +
  scale_color_manual(values = c("Prior" = "red", "Posterior" = "blue")) +
  scale_fill_manual(values = c("Prior" = "red", "Posterior" = "blue")) +
  theme_minimal() +
  theme(legend.title = element_blank())

# Plot Improvement Ratio
p2 <- ggplot(results_df, aes(x = prop, y = improvement.mean)) +
  geom_line(color = "green", linewidth = 1) +
  geom_point(size = 2, color = "green") +
  geom_ribbon(aes(ymin = improvement.mean - improvement.sd, 
                  ymax = improvement.mean + improvement.sd), 
              alpha = 0.2, fill = "green") +
  scale_x_continuous(labels = percent, limits = c(xlim_min, 1)) +
  labs(title = "Improvement in Calibration",
       x = "Proportion of Shared Populations",
       y = "Relative Improvement in the") +
  theme_minimal()

# Ensure the directory exists
if (!dir.exists("figures")) {
  dir.create("figures")
}

# Save plots
ggsave(paste0("figures/the_prior_vs_posterior_", data_source, ".jpg"),
       p1, width = 10, height = 6, dpi = 300)
ggsave(paste0("figures/improvement_ratio_", data_source, ".jpg"),
       p2, width = 10, height = 6, dpi = 300)

cat("Visualization completed. Plots saved.\n")
