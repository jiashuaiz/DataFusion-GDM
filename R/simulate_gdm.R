# Globals for NSE in ggplot2
utils::globalVariables(c("Dim1", "Dim2", "MajorGroup", "Bottlenecked", "Admixed", "Col", "Row", "Distance"))

# Simulation and visualization of genetic distances

# Simulate genetic distances using realistic population structure
# (roxygen documentation is placed above simulate_genetic_distances())
# Utility functions
create_population_metadata <- function(
    n_pops = 50,
    n_major_groups = 5, 
    n_subgroups = 12,
    use_subgroups = TRUE
) {
  group_probs <- stats::rbeta(n_major_groups, 2, 2)
  group_probs <- group_probs / sum(group_probs)
  pop_groups <- sample(1:n_major_groups, n_pops, replace = TRUE, prob = group_probs)
  pop_labels <- paste0("Pop", sprintf("%02d", seq_len(n_pops)))
  pop_info <- data.frame(
    Population = pop_labels,
    MajorGroup = paste0("Group", pop_groups)
  )
  if (use_subgroups) {
    subgroup_assignment <- rep(0, n_pops)
    subgroups_per_group <- ceiling(n_subgroups / n_major_groups)
    for (g in seq_len(n_major_groups)) {
      group_pops <- which(pop_groups == g)
      if (length(group_pops) > 0) {
        start_subgroup <- (g-1) * subgroups_per_group + 1
        end_subgroup <- min(g * subgroups_per_group, n_subgroups)
        available_subgroups <- start_subgroup:end_subgroup
        if (length(available_subgroups) > 1) {
          subgroup_probs <- stats::rbeta(length(available_subgroups), 1.5, 1.5)
          subgroup_probs <- subgroup_probs / sum(subgroup_probs)
          subgroup_assignment[group_pops] <- sample(
            available_subgroups, 
            length(group_pops), 
            replace = TRUE, 
            prob = subgroup_probs
          )
        } else {
          subgroup_assignment[group_pops] <- available_subgroups
        }
      }
    }
    pop_info$Subgroup <- paste0("Subgroup", subgroup_assignment)
  } else {
    pop_info$Subgroup <- paste0("Subgroup", seq_len(n_pops))
  }
  return(pop_info)
}

generate_group_centroids <- function(
    n_major_groups,
    n_dims,
    geo_dims,
    group_separation = 15,
    use_genetic_dims = TRUE
) {
  group_centroids <- matrix(0, nrow = n_major_groups, ncol = n_dims)
  for (g in seq_len(n_major_groups)) {
    if (geo_dims > 0) {
      if (g <= geo_dims + 1) {
        angle <- 2 * pi * (g-1) / (geo_dims + 1)
        group_separation_vec <- c(cos(angle), sin(angle))
        group_centroids[g, 1:geo_dims] <- group_separation * group_separation_vec[1:geo_dims]
      } else {
        group_centroids[g, 1:geo_dims] <- group_separation * stats::runif(geo_dims, -1, 1)
      }
    }
    if (use_genetic_dims && n_dims > geo_dims) {
      group_centroids[g, (geo_dims+1):n_dims] <- group_separation * 0.5 * stats::rnorm(n_dims - geo_dims)
    }
  }
  return(group_centroids)
}

generate_subgroup_centroids <- function(
    n_subgroups,
    n_major_groups,
    group_centroids,
    n_dims,
    subgroup_separation = 5,
    use_subgroups = TRUE
) {
  if (!use_subgroups) {
    return(group_centroids)
  }
  subgroup_centroids <- matrix(0, nrow = n_subgroups, ncol = n_dims)
  subgroups_per_group <- ceiling(n_subgroups / n_major_groups)
  for (s in seq_len(n_subgroups)) {
    parent_group <- ceiling(s / subgroups_per_group)
    if (parent_group > n_major_groups) parent_group <- n_major_groups
    subgroup_centroids[s,] <- group_centroids[parent_group,] + 
      subgroup_separation * stats::rnorm(n_dims, 0, 0.5)
  }
  return(subgroup_centroids)
}

position_populations <- function(
    n_pops,
    pop_info,
    subgroup_centroids,
    n_dims,
    pop_dispersion = 0.5
) {
  pop_positions <- matrix(0, nrow = n_pops, ncol = n_dims)
  for (i in seq_len(n_pops)) {
    sg_str <- pop_info$Subgroup[i]
    sg <- as.numeric(gsub("Subgroup", "", sg_str))
    if (sg > nrow(subgroup_centroids) || sg < 1) {
      sg <- 1
    }
    pop_positions[i,] <- subgroup_centroids[sg,] + pop_dispersion * stats::rnorm(n_dims, 0, 1)
  }
  return(pop_positions)
}

create_admixed_populations <- function(
    pop_positions,
    pop_info,
    n_pops,
    subgroup_centroids,
    n_dims,
    admixture_prob = 0.1,
    pop_dispersion = 0.5,
    use_admixture = TRUE
) {
  pop_info$Admixed <- "No"
  if (!use_admixture) {
    return(list(positions = pop_positions, metadata = pop_info))
  }
  n_admixed <- round(n_pops * admixture_prob)
  if (n_admixed > 0 && nrow(subgroup_centroids) > 1) {
    admixed_pops <- sample(seq_len(n_pops), n_admixed)
    for (i in admixed_pops) {
      sg_str <- pop_info$Subgroup[i]
      sg1 <- as.numeric(gsub("Subgroup", "", sg_str))
      if (sg1 > nrow(subgroup_centroids) || sg1 < 1) {
        sg1 <- 1
      }
      all_sgs <- seq_len(nrow(subgroup_centroids))
      if (length(all_sgs) > 1) {
        sg2 <- sample(all_sgs[all_sgs != sg1], 1)
        admix_ratio <- stats::runif(1, 0.3, 0.7)
        pop_positions[i,] <- admix_ratio * subgroup_centroids[sg1,] + 
          (1-admix_ratio) * subgroup_centroids[sg2,] +
          pop_dispersion * 0.5 * stats::rnorm(n_dims, 0, 1)
        pop_info$Admixed[i] <- "Yes"
      }
    }
  }
  return(list(positions = pop_positions, metadata = pop_info))
}

create_bottlenecked_populations <- function(
    pop_positions,
    pop_info,
    n_pops,
    subgroup_centroids,
    geo_dims,
    n_dims,
    bottleneck_prob = 0.05,
    pop_dispersion = 0.5,
    use_bottlenecks = TRUE
) {
  pop_info$Bottlenecked <- "No"
  if (!use_bottlenecks) {
    return(list(positions = pop_positions, metadata = pop_info))
  }
  n_bottleneck <- round(n_pops * bottleneck_prob)
  if (n_bottleneck > 0) {
    admixed_pops <- which(pop_info$Admixed == "Yes")
    candidate_pops <- setdiff(seq_len(n_pops), admixed_pops)
    if (length(candidate_pops) < n_bottleneck) {
      candidate_pops <- seq_len(n_pops)
    }
    bottleneck_pops <- sample(candidate_pops, min(n_bottleneck, length(candidate_pops)))
    for (i in bottleneck_pops) {
      sg_str <- pop_info$Subgroup[i]
      sg <- as.numeric(gsub("Subgroup", "", sg_str))
      if (sg > nrow(subgroup_centroids) || sg < 1) {
        sg <- 1
      }
      pop_positions[i,] <- subgroup_centroids[sg,] + 0.3 * pop_dispersion * stats::rnorm(n_dims, 0, 1)
      if (n_dims > geo_dims) {
        genetic_dims_count <- n_dims - geo_dims
        pop_positions[i, (geo_dims+1):n_dims] <- pop_positions[i, (geo_dims+1):n_dims] + 
          stats::rnorm(genetic_dims_count, 0, 1)
      }
      pop_info$Bottlenecked[i] <- "Yes"
    }
  }
  return(list(positions = pop_positions, metadata = pop_info))
}

calculate_raw_distances <- function(
    pop_positions,
    n_pops,
    geo_dims,
    isolation_factor = 0.8,
    use_isolation_by_distance = TRUE
) {
  dist_matrix <- matrix(0, nrow = n_pops, ncol = n_pops)
  for (i in seq_len(n_pops)) {
    for (j in i:n_pops) {
      if (i != j) {
        geo_dist <- if (geo_dims > 0) sqrt(sum((pop_positions[i, 1:geo_dims] - pop_positions[j, 1:geo_dims])^2)) else 0
        gen_dist <- 0
        if (ncol(pop_positions) > geo_dims) {
          gen_dims <- (geo_dims+1):ncol(pop_positions)
          gen_dist <- sqrt(sum((pop_positions[i, gen_dims] - pop_positions[j, gen_dims])^2))
        }
        if (use_isolation_by_distance) {
          raw_dist <- isolation_factor * geo_dist + (1 - isolation_factor) * gen_dist
        } else {
          raw_dist <- sqrt(sum((pop_positions[i,] - pop_positions[j,])^2))
        }
        dist_matrix[i, j] <- raw_dist
        dist_matrix[j, i] <- raw_dist
      }
    }
  }
  return(dist_matrix)
}

check_distance_matrix <- function(dist_matrix, message = "Distance matrix check", verbose = TRUE) {
  if (!verbose) return()
  cat("\n", message, ":\n")
  cat("Dimensions:", dim(dist_matrix), "\n")
  cat("Is symmetric:", all(abs(dist_matrix - t(dist_matrix)) < 1e-10), "\n")
  rng <- range(dist_matrix)
  cat("Range of distances:", rng, "\n")
  cat("Mean distance:", mean(dist_matrix), "\n")
  if (stats::sd(as.vector(dist_matrix)) < 1e-10) {
    cat("WARNING: All distances are essentially identical!\n")
  }
  diag_vals <- diag(dist_matrix)
  cat("Diagonal values - min:", min(diag_vals), "max:", max(diag_vals), "\n")
  quants <- stats::quantile(as.vector(dist_matrix))
  cat("Quantiles: 0%:", quants[1], "25%:", quants[2], "50%:", quants[3], 
      "75%:", quants[4], "100%:", quants[5], "\n")
}

transform_to_genetic_distances <- function(
    raw_distances,
    n_pops,
    nonlinear_factor = 0.7,
    noise_level = 0.1,
    use_nonlinear = TRUE,
    use_noise = TRUE,
    verbose = TRUE
) {
  check_distance_matrix(raw_distances, "Raw distances before transformation", verbose)
  genetic_dist <- matrix(0, nrow = n_pops, ncol = n_pops)
  for (i in seq_len(n_pops)) {
    for (j in i:n_pops) {
      if (i != j) {
        raw_dist <- raw_distances[i, j]
        if (use_nonlinear) {
          dist_val <- raw_dist^nonlinear_factor / (1 + raw_dist^nonlinear_factor)
        } else {
          dist_val <- raw_dist / max(raw_distances)
        }
        if (use_noise) {
          noise_factor <- noise_level * (0.5 + dist_val)
          dist_val <- dist_val * (1 + stats::rnorm(1, 0, noise_factor))
        }
        dist_val <- max(0, min(1, dist_val))
        genetic_dist[i, j] <- dist_val
        genetic_dist[j, i] <- dist_val
      }
    }
  }
  check_distance_matrix(genetic_dist, "Final genetic distances after transformation", verbose)
  return(genetic_dist)
}

#' Create a heatmap of genetic distances (ggplot2)
#' @description Returns a ggplot heatmap of the distance matrix using ggplot2 only (no Bioconductor dependencies).
#' @param dist_matrix Symmetric numeric distance matrix with row/column names
#' @param pop_info Data frame with at least `Population` and `MajorGroup` columns
#' @return A ggplot object
create_distance_heatmap <- function(dist_matrix, pop_info) {
  rn <- rownames(dist_matrix); cn <- colnames(dist_matrix)
  if (is.null(rn) || is.null(cn)) {
    stop("Distance matrix must have row and column names")
  }
  df <- as.data.frame(as.table(dist_matrix))
  names(df) <- c("Row", "Col", "Distance")
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Col, y = Row, fill = Distance)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "white", high = "#1f77b4") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(x = NULL, y = NULL)
  p
}

#' Create MDS plot of genetic distances
#' @param dist_matrix Symmetric numeric distance matrix
#' @param pop_info Data frame with metadata columns
#' @return A ggplot object
create_mds_plot <- function(dist_matrix, pop_info) {
  rownames(dist_matrix) <- pop_info$Population
  colnames(dist_matrix) <- pop_info$Population
  dist_obj <- stats::as.dist(dist_matrix)
  mds_coords <- as.data.frame(stats::cmdscale(dist_obj, k = 2))
  colnames(mds_coords) <- c("Dim1", "Dim2")
  mds_coords$Population <- rownames(mds_coords)
  mds_coords$MajorGroup <- pop_info$MajorGroup[match(mds_coords$Population, pop_info$Population)]
  mds_coords$Subgroup <- pop_info$Subgroup[match(mds_coords$Population, pop_info$Population)]
  mds_coords$Admixed <- pop_info$Admixed[match(mds_coords$Population, pop_info$Population)]
  mds_coords$Bottlenecked <- pop_info$Bottlenecked[match(mds_coords$Population, pop_info$Population)]
  p <- ggplot2::ggplot(mds_coords, ggplot2::aes(x = Dim1, y = Dim2, color = MajorGroup)) +
    ggplot2::geom_point(ggplot2::aes(shape = ifelse(Bottlenecked == "Yes", "Bottlenecked", "Normal")), 
               size = 3, alpha = 0.8) +
    ggplot2::scale_shape_manual(values = c("Bottlenecked" = 17, "Normal" = 16)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "MDS Plot of Simulated Genetic Distances",
         x = "Dimension 1", y = "Dimension 2",
         shape = "Population Type") +
    ggplot2::theme(legend.position = "right")
  if (any(pop_info$Admixed == "Yes")) {
    p <- p + ggplot2::geom_point(data = subset(mds_coords, Admixed == "Yes"), 
                 size = 5, shape = 1, color = "black")
  }
  return(p)
}

#' Simulate genetic distances using realistic population structure
#' @description Generates a synthetic genetic distance matrix and metadata using hierarchical population structure, admixture and bottleneck options.
#' @param n_pops Number of populations
#' @param n_major_groups Number of major groups
#' @param n_subgroups Number of subgroups
#' @param geo_dims Geographic dimensions
#' @param genetic_dims Additional genetic drift dimensions
#' @param group_separation Separation between major groups
#' @param subgroup_separation Separation between subgroups
#' @param pop_dispersion Within-subgroup dispersion
#' @param isolation_factor Weight for geography in isolation-by-distance model (0-1)
#' @param admixture_prob Proportion of admixed populations
#' @param bottleneck_prob Proportion of bottlenecked populations
#' @param noise_level Noise level in transformation
#' @param nonlinear_factor Nonlinearity factor in transformation
#' @param use_subgroups Whether to create subgroups
#' @param use_genetic_dims Whether to include genetic dimensions
#' @param use_admixture Whether to include admixture
#' @param use_bottlenecks Whether to include bottlenecks
#' @param use_isolation_by_distance Whether to weight geographic distance
#' @param use_nonlinear Whether to apply nonlinear transformation
#' @param use_noise Whether to add noise
#' @param seed Optional seed for reproducibility (NULL leaves the RNG state unchanged)
#' @param verbose Print diagnostics
#' @return A list with `distance_matrix`, `population_info`, `position_matrix`, and `parameters`.
#' @export
simulate_genetic_distances <- function(
    n_pops = 50,
    n_major_groups = 5,
    n_subgroups = 12,
    geo_dims = 2,
    genetic_dims = 2,
    group_separation = 15,
    subgroup_separation = 5,
    pop_dispersion = 0.5,
    isolation_factor = 0.8,
    admixture_prob = 0.1,
    bottleneck_prob = 0.05,
    noise_level = 0.1,
    nonlinear_factor = 0.7,
    use_subgroups = TRUE,
    use_genetic_dims = TRUE,
    use_admixture = TRUE,
    use_bottlenecks = TRUE,
    use_isolation_by_distance = TRUE,
    use_nonlinear = TRUE,
  use_noise = TRUE,
  seed = NULL,
    verbose = TRUE
) {
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed)) {
      stop("`seed` must be a single finite numeric value or NULL.")
    }
    seed_backup_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (seed_backup_exists) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    on.exit({
      if (seed_backup_exists) {
        if (is.null(old_seed)) {
          if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
          }
        } else {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        }
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(seed))
  }
  n_pops <- max(10, n_pops)
  n_major_groups <- min(n_pops/2, max(2, n_major_groups))
  n_subgroups <- min(n_pops, max(n_major_groups, n_subgroups))
  if (!use_genetic_dims) genetic_dims <- 0
  n_dims <- geo_dims + genetic_dims
  pop_info <- create_population_metadata(n_pops, n_major_groups, n_subgroups, use_subgroups)
  group_centroids <- generate_group_centroids(n_major_groups, n_dims, geo_dims, group_separation, use_genetic_dims)
  subgroup_centroids <- generate_subgroup_centroids(n_subgroups, n_major_groups, group_centroids, n_dims, subgroup_separation, use_subgroups)
  pop_positions <- position_populations(n_pops, pop_info, subgroup_centroids, n_dims, pop_dispersion)
  admix_result <- create_admixed_populations(pop_positions, pop_info, n_pops, subgroup_centroids, n_dims, admixture_prob, pop_dispersion, use_admixture)
  pop_positions <- admix_result$positions
  pop_info <- admix_result$metadata
  bottleneck_result <- create_bottlenecked_populations(pop_positions, pop_info, n_pops, subgroup_centroids, geo_dims, n_dims, bottleneck_prob, pop_dispersion, use_bottlenecks)
  pop_positions <- bottleneck_result$positions
  pop_info <- bottleneck_result$metadata
  raw_distances <- calculate_raw_distances(pop_positions, n_pops, geo_dims, isolation_factor, use_isolation_by_distance)
  genetic_distances <- transform_to_genetic_distances(raw_distances, n_pops, nonlinear_factor, noise_level, use_nonlinear, use_noise, verbose)
  genetic_distances <- (genetic_distances + t(genetic_distances)) / 2
  rownames(genetic_distances) <- pop_info$Population
  colnames(genetic_distances) <- pop_info$Population
  list(
    distance_matrix = genetic_distances,
    population_info = pop_info,
    position_matrix = pop_positions,
    parameters = list(
      n_pops = n_pops,
      n_major_groups = n_major_groups,
      n_subgroups = n_subgroups,
      geo_dims = geo_dims,
      genetic_dims = genetic_dims,
      group_separation = group_separation,
      subgroup_separation = subgroup_separation,
      pop_dispersion = pop_dispersion,
      isolation_factor = isolation_factor,
      admixture_prob = admixture_prob,
      bottleneck_prob = bottleneck_prob,
      use_subgroups = use_subgroups,
      use_genetic_dims = use_genetic_dims,
      use_admixture = use_admixture,
      use_bottlenecks = use_bottlenecks,
      use_isolation_by_distance = use_isolation_by_distance,
      use_nonlinear = use_nonlinear,
      use_noise = use_noise,
      seed = seed
    )
  )
}

#' Create plotting handles for simulation results
#' @param sim_results A list returned by simulate_genetic_distances() or run_genetic_simulation()
#' @return A list with `heatmap` and `mds` functions that print plots when called
#' @export
visualize_results <- function(sim_results) {
  dist_matrix <- sim_results$distance_matrix
  pop_info <- sim_results$population_info
  # Always create MDS plot
  mds_plot <- create_mds_plot(dist_matrix, pop_info)
  # Simple heatmap via ggplot2 to avoid non-CRAN Suggests
  heatmap_obj <- create_distance_heatmap(dist_matrix, pop_info)
  list(
    heatmap = function() { print(heatmap_obj) },
    mds = function() { print(mds_plot) }
  )
}

#' Run a high-level genetic simulation with configurable model
#' @inheritParams simulate_genetic_distances
#' @param model One of "mixed", "geography", "genetics", or "custom"
#' @param geo_dims Geographic dimensions (overrides default based on model if set)
#' @param isolation_factor Geography-genetics balance (overrides default based on model if set)
#' @param genetic_dims Genetic dimensions (overrides default based on model if set)
#' @param subgroup_separation Separation between subgroups (default: group_separation/3 when NULL)
#' @param output_file Optional CSV file path to write the distance matrix
#' @param seed Optional seed forwarded to simulate_genetic_distances()
#' @return List with `results` and `plots` (functions to print plots)
#' @export
run_genetic_simulation <- function(
  n_pops = 30,
  n_major_groups = 4,
  n_subgroups = 8,
  model = "mixed",
  geo_dims = NULL,
  isolation_factor = NULL,
  genetic_dims = NULL,
  group_separation = 15,
  subgroup_separation = NULL,
  pop_dispersion = 0.5,
  admixture_prob = 0.15,
  bottleneck_prob = 0.10,
  use_subgroups = TRUE,
  use_genetic_dims = NULL,
  use_admixture = TRUE,
  use_bottlenecks = TRUE,
  use_isolation_by_distance = NULL,
  use_nonlinear = TRUE,
  use_noise = TRUE,
  seed = NULL,
  output_file = NULL,
  verbose = TRUE
) {
  if (model == "geography") {
    geo_dims <- ifelse(is.null(geo_dims), 2, geo_dims)
    genetic_dims <- 0
    isolation_factor <- 1.0
    use_isolation_by_distance <- TRUE
    use_genetic_dims <- FALSE
  } else if (model == "genetics") {
    geo_dims <- 0
    genetic_dims <- ifelse(is.null(genetic_dims), 2, genetic_dims)
    isolation_factor <- 0.0
    use_isolation_by_distance <- FALSE
    use_genetic_dims <- TRUE
  } else if (model == "mixed") {
    geo_dims <- ifelse(is.null(geo_dims), 2, geo_dims)
    genetic_dims <- ifelse(is.null(genetic_dims), 2, genetic_dims)
    isolation_factor <- ifelse(is.null(isolation_factor), 0.8, isolation_factor)
    use_isolation_by_distance <- TRUE
    use_genetic_dims <- TRUE
  } else if (model == "custom") {
    geo_dims <- ifelse(is.null(geo_dims), 2, geo_dims)
    genetic_dims <- ifelse(is.null(genetic_dims), 2, genetic_dims)
    isolation_factor <- ifelse(is.null(isolation_factor), 0.8, isolation_factor)
    use_isolation_by_distance <- ifelse(is.null(use_isolation_by_distance), isolation_factor > 0, use_isolation_by_distance)
    use_genetic_dims <- ifelse(is.null(use_genetic_dims), genetic_dims > 0, use_genetic_dims)
  } else {
    stop("Invalid model specified. Choose from: 'mixed', 'geography', 'genetics', or 'custom')")
  }
  if (use_isolation_by_distance && is.null(isolation_factor)) isolation_factor <- 0.8
  if (isolation_factor == 1.0) use_isolation_by_distance <- TRUE
  if (geo_dims == 0) { use_isolation_by_distance <- FALSE; isolation_factor <- 0 }
  if (genetic_dims == 0) { use_genetic_dims <- FALSE }
  if (use_genetic_dims && genetic_dims == 0) { genetic_dims <- 2 }
  if (is.null(subgroup_separation)) subgroup_separation <- group_separation / 3
  sim_results <- simulate_genetic_distances(
    n_pops = n_pops,
    n_major_groups = n_major_groups,
    n_subgroups = n_subgroups,
    geo_dims = geo_dims,
    genetic_dims = genetic_dims,
    group_separation = group_separation,
    subgroup_separation = subgroup_separation,
    pop_dispersion = pop_dispersion,
    isolation_factor = isolation_factor,
    admixture_prob = admixture_prob,
    bottleneck_prob = bottleneck_prob,
    use_subgroups = use_subgroups,
    use_genetic_dims = use_genetic_dims,
    use_admixture = use_admixture,
    use_bottlenecks = use_bottlenecks,
    use_isolation_by_distance = use_isolation_by_distance,
    use_nonlinear = use_nonlinear,
    use_noise = use_noise,
    seed = seed,
    verbose = verbose
  )
  plots <- visualize_results(sim_results)
  if (!is.null(output_file)) {
    csv_output <- if (grepl("\\.csv$", output_file)) output_file else paste0(output_file, ".csv")
    utils::write.csv(sim_results$distance_matrix, csv_output)
  }
  list(results = sim_results, plots = plots)
}

#' Run simulation with predefined biological scenarios
#' @param scenario Scenario name: 'default', 'island', 'stepping_stone', 'admixture', 'ancient_divergence', 'simple'
#' @param n_pops Number of populations
#' @param output_file Optional CSV path to write the distance matrix
#' @param seed Optional seed forwarded to run_genetic_simulation()
#' @param verbose Print diagnostic information
#' @return Same structure as run_genetic_simulation()
#' @export
run_genetic_scenario <- function(
    scenario = "default",
    n_pops = 30,
    output_file = NULL,
    seed = NULL,
    verbose = TRUE
) {
  scenarios <- list(
    default = list(
      model = "mixed", n_major_groups = 4, n_subgroups = 8, admixture_prob = 0.15, bottleneck_prob = 0.1,
      group_separation = 15, subgroup_separation = 5, pop_dispersion = 0.5,
      use_subgroups = TRUE, use_admixture = TRUE, use_bottlenecks = TRUE, use_nonlinear = TRUE, use_noise = TRUE
    ),
    island = list(
      model = "mixed", n_major_groups = 6, n_subgroups = 6, admixture_prob = 0.05, bottleneck_prob = 0.2,
      group_separation = 25, subgroup_separation = 0, pop_dispersion = 0.3,
      use_subgroups = FALSE, use_admixture = TRUE, use_bottlenecks = TRUE, use_nonlinear = TRUE, use_noise = TRUE
    ),
    stepping_stone = list(
      model = "geography", n_major_groups = 3, n_subgroups = 12, admixture_prob = 0.1, bottleneck_prob = 0.05,
      group_separation = 10, subgroup_separation = 3, pop_dispersion = 0.7,
      use_subgroups = TRUE, use_admixture = TRUE, use_bottlenecks = FALSE, use_nonlinear = TRUE, use_noise = TRUE
    ),
    admixture = list(
      model = "mixed", n_major_groups = 3, n_subgroups = 6, admixture_prob = 0.4, bottleneck_prob = 0,
      group_separation = 12, subgroup_separation = 4, pop_dispersion = 0.6,
      use_subgroups = TRUE, use_admixture = TRUE, use_bottlenecks = FALSE, use_nonlinear = TRUE, use_noise = TRUE
    ),
    ancient_divergence = list(
      model = "genetics", n_major_groups = 4, n_subgroups = 10, admixture_prob = 0.05, bottleneck_prob = 0.1,
      group_separation = 30, subgroup_separation = 8, pop_dispersion = 0.4,
      use_subgroups = TRUE, use_admixture = TRUE, use_bottlenecks = TRUE, use_nonlinear = TRUE, use_noise = FALSE
    ),
    simple = list(
      model = "custom", geo_dims = 2, genetic_dims = 0, n_major_groups = 3, n_subgroups = 0,
      admixture_prob = 0, bottleneck_prob = 0, group_separation = 15, subgroup_separation = 0, pop_dispersion = 0.5,
      use_subgroups = FALSE, use_admixture = FALSE, use_bottlenecks = FALSE, use_nonlinear = FALSE, use_noise = FALSE
    )
  )
  if (!scenario %in% names(scenarios)) stop("Unknown scenario: ", scenario)
  params <- scenarios[[scenario]]
  if (!is.null(output_file)) {
    output_file <- paste0(output_file, "_", scenario)
  }
  do.call(run_genetic_simulation, c(list(n_pops = n_pops, output_file = output_file, seed = seed, verbose = verbose), params))
}

#' Export a simulated GDM to CSV
#' @param output_file Output CSV filename (defaults to a session-scoped temporary path)
#' @param scenario Scenario name
#' @param n_pops Number of populations
#' @param verbose Verbose output
#' @param seed Optional seed forwarded to run_genetic_scenario()
#' @return Invisibly, the normalized path to the written CSV
#' @export
#' @examples
#' tmp <- export_simulated_gdm(verbose = FALSE)
#' if (file.exists(tmp)) unlink(tmp)
export_simulated_gdm <- function(output_file = tempfile("gdm_", fileext = ".csv"),
                                 scenario = "default",
                                 n_pops = 30,
                                 verbose = TRUE,
                                 seed = NULL) {
  if (is.null(output_file) || !nzchar(output_file)) {
    stop("`output_file` must be a non-empty string path.")
  }
  target_dir <- dirname(output_file)
  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  }
  # Normalize without extension for internal write (run_genetic_scenario appends _scenario)
  base_no_ext <- sub("\\.csv$", "", output_file)
  # Run and write
  invisible(run_genetic_scenario(scenario = scenario, n_pops = n_pops, output_file = base_no_ext, seed = seed, verbose = verbose))
  # Compute actual output path used by run_genetic_simulation
  actual <- paste0(base_no_ext, "_", scenario, ".csv")
  # Return normalized existing path if possible
  out <- tryCatch(normalizePath(actual, mustWork = TRUE), error = function(e) actual)
  invisible(out)
}
