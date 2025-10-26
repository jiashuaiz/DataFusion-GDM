####################################################
# GENETIC DISTANCE SIMULATION
####################################################

# 1. UTILITY FUNCTIONS ####################################################

#' Create population metadata
#' return Data frame with population group assignments
create_population_metadata <- function(
    n_pops = 50,
    n_major_groups = 5, 
    n_subgroups = 12,
    use_subgroups = TRUE
) {
  # Create uneven group sizes for realism
  group_probs <- rbeta(n_major_groups, 2, 2)
  group_probs <- group_probs / sum(group_probs)
  
  # Assign populations to major groups
  pop_groups <- sample(1:n_major_groups, n_pops, replace = TRUE, prob = group_probs)
  
  # Create population labels
  pop_labels <- paste0("Pop", sprintf("%02d", 1:n_pops))
  
  # Initialize metadata
  pop_info <- data.frame(
    Population = pop_labels,
    MajorGroup = paste0("Group", pop_groups)
  )
  
  # Assign subgroups only if requested
  if (use_subgroups) {
    subgroup_assignment <- rep(0, n_pops)
    subgroups_per_group <- ceiling(n_subgroups / n_major_groups)
    
    for (g in 1:n_major_groups) {
      group_pops <- which(pop_groups == g)
      if (length(group_pops) > 0) {
        # Create subgroups within this major group
        start_subgroup <- (g-1) * subgroups_per_group + 1
        end_subgroup <- min(g * subgroups_per_group, n_subgroups)
        available_subgroups <- start_subgroup:end_subgroup
        
        # Check if there's more than one subgroup before using probability weights
        if (length(available_subgroups) > 1) {
          subgroup_probs <- rbeta(length(available_subgroups), 1.5, 1.5)
          subgroup_probs <- subgroup_probs / sum(subgroup_probs)
          
          subgroup_assignment[group_pops] <- sample(
            available_subgroups, 
            length(group_pops), 
            replace = TRUE, 
            prob = subgroup_probs
          )
        } else {
          # If only one subgroup is available, assign all populations to it
          subgroup_assignment[group_pops] <- available_subgroups
        }
      }
    }
    
    pop_info$Subgroup <- paste0("Subgroup", subgroup_assignment)
  } else {
    # If not using subgroups, each population is its own subgroup
    pop_info$Subgroup <- paste0("Subgroup", 1:n_pops)
  }
  
  return(pop_info)
}

#' Generate group centroids in multidimensional space
#' return Matrix of group centroid coordinates
generate_group_centroids <- function(
    n_major_groups,
    n_dims,
    geo_dims,
    group_separation = 15,
    use_genetic_dims = TRUE
) {
  group_centroids <- matrix(0, nrow = n_major_groups, ncol = n_dims)
  
  # Generate major group centroids that are well-separated
  for (g in 1:n_major_groups) {
    # Only place geographic centroids if geo_dims > 0
    if (geo_dims > 0) {
      # Place centroids roughly in a circle for the geographic dimensions
      if (g <= geo_dims + 1) {
        angle <- 2 * pi * (g-1) / (geo_dims + 1)
        group_centroids[g, 1:geo_dims] <- group_separation * c(cos(angle), sin(angle))[1:geo_dims]
      } else {
        # Random placement for additional groups
        group_centroids[g, 1:geo_dims] <- group_separation * runif(geo_dims, -1, 1)
      }
    }
    
    # For genetic dimensions, add random values if requested
    if (use_genetic_dims && n_dims > geo_dims) {
      group_centroids[g, (geo_dims+1):n_dims] <- group_separation * 0.5 * rnorm(n_dims - geo_dims)
    }
  }
  
  return(group_centroids)
}

#' Generate subgroup centroids around major group centroids
#' return Matrix of subgroup centroid coordinates
generate_subgroup_centroids <- function(
    n_subgroups,
    n_major_groups,
    group_centroids,
    n_dims,
    subgroup_separation = 5,
    use_subgroups = TRUE
) {
  if (!use_subgroups) {
    # If not using subgroups, return the major group centroids
    return(group_centroids)
  }
  
  subgroup_centroids <- matrix(0, nrow = n_subgroups, ncol = n_dims)
  subgroups_per_group <- ceiling(n_subgroups / n_major_groups)
  
  for (s in 1:n_subgroups) {
    # Determine which major group this subgroup belongs to
    parent_group <- ceiling(s / subgroups_per_group)
    if (parent_group > n_major_groups) parent_group <- n_major_groups
    
    # Place subgroup near its parent group centroid
    subgroup_centroids[s,] <- group_centroids[parent_group,] + 
      subgroup_separation * rnorm(n_dims, 0, 0.5)
  }
  
  return(subgroup_centroids)
}

# 2. POPULATION POSITION FUNCTIONS ####################################################

#' Position populations in multidimensional space
#' return Matrix of population coordinates
position_populations <- function(
    n_pops,
    pop_info,
    subgroup_centroids,
    n_dims,
    pop_dispersion = 0.5
) {
  # Initialize position matrix
  pop_positions <- matrix(0, nrow = n_pops, ncol = n_dims)
  
  # Position each population near its subgroup centroid
  for (i in 1:n_pops) {
    # Extract the subgroup ID correctly
    sg_str <- pop_info$Subgroup[i]
    sg <- as.numeric(gsub("Subgroup", "", sg_str))
    
    # Handle case when subgroup number is out of bounds
    if (sg > nrow(subgroup_centroids) || sg < 1) {
      sg <- 1
    }
    
    pop_positions[i,] <- subgroup_centroids[sg,] + pop_dispersion * rnorm(n_dims, 0, 1)
  }
  
  return(pop_positions)
}

#' Create admixed populations
#' return Updated population positions and metadata
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
  # Create a column for admixture status
  pop_info$Admixed <- "No"
  
  # Skip if admixture is disabled
  if (!use_admixture) {
    return(list(
      positions = pop_positions,
      metadata = pop_info
    ))
  }
  
  # Identify populations for admixture
  n_admixed <- round(n_pops * admixture_prob)
  
  if (n_admixed > 0 && nrow(subgroup_centroids) > 1) {
    admixed_pops <- sample(1:n_pops, n_admixed)
    
    for (i in admixed_pops) {
      # Extract the subgroup ID correctly
      sg_str <- pop_info$Subgroup[i]
      sg1 <- as.numeric(gsub("Subgroup", "", sg_str))
      
      # Handle case when subgroup number is out of bounds
      if (sg1 > nrow(subgroup_centroids) || sg1 < 1) {
        sg1 <- 1
      }
      
      # Choose a different subgroup for admixture
      all_sgs <- 1:nrow(subgroup_centroids)
      if (length(all_sgs) > 1) {
        sg2 <- sample(all_sgs[all_sgs != sg1], 1)
        
        # Create admixed position - weighted average of the two subgroup positions
        admix_ratio <- runif(1, 0.3, 0.7)  # random admixture proportion
        pop_positions[i,] <- admix_ratio * subgroup_centroids[sg1,] + 
          (1-admix_ratio) * subgroup_centroids[sg2,] +
          pop_dispersion * 0.5 * rnorm(n_dims, 0, 1)  # reduced noise for admixed pops
        
        # Update metadata
        pop_info$Admixed[i] <- "Yes"
      }
    }
  }
  
  return(list(
    positions = pop_positions,
    metadata = pop_info
  ))
}

#' Create bottlenecked populations
#' return Updated population positions and metadata
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
  # Create a column for bottleneck status
  pop_info$Bottlenecked <- "No"
  
  # Skip if bottlenecks are disabled
  if (!use_bottlenecks) {
    return(list(
      positions = pop_positions,
      metadata = pop_info
    ))
  }
  
  # Identify populations for bottleneck effects
  n_bottleneck <- round(n_pops * bottleneck_prob)
  
  if (n_bottleneck > 0) {
    # Exclude already admixed populations if possible
    admixed_pops <- which(pop_info$Admixed == "Yes")
    candidate_pops <- setdiff(1:n_pops, admixed_pops)
    
    # If all populations are admixed, we'll have to use some of them
    if (length(candidate_pops) < n_bottleneck) {
      candidate_pops <- 1:n_pops
    }
    
    bottleneck_pops <- sample(candidate_pops, min(n_bottleneck, length(candidate_pops)))
    
    for (i in bottleneck_pops) {
      # Extract the subgroup ID correctly
      sg_str <- pop_info$Subgroup[i]
      sg <- as.numeric(gsub("Subgroup", "", sg_str))
      
      # Handle case when subgroup number is out of bounds
      if (sg > nrow(subgroup_centroids) || sg < 1) {
        sg <- 1
      }
      
      # Bottlenecked populations have reduced genetic diversity (closer to subgroup centroid)
      # but might be slightly offset due to drift
      pop_positions[i,] <- subgroup_centroids[sg,] + 
        0.3 * pop_dispersion * rnorm(n_dims, 0, 1)  # reduced dispersion
      
      # Add additional drift in genetic dimensions if they exist
      if (n_dims > geo_dims) {
        genetic_dims_count <- n_dims - geo_dims
        pop_positions[i, (geo_dims+1):n_dims] <- pop_positions[i, (geo_dims+1):n_dims] + 
          rnorm(genetic_dims_count, 0, 1)
      }
      
      # Update metadata
      pop_info$Bottlenecked[i] <- "Yes"
    }
  }
  
  return(list(
    positions = pop_positions,
    metadata = pop_info
  ))
}

# 3. DISTANCE CALCULATION FUNCTIONS ####################################################

#' Calculate raw Euclidean distances between populations
#' return Basic distance matrix
calculate_raw_distances <- function(
    pop_positions,
    n_pops,
    geo_dims,
    isolation_factor = 0.8,
    use_isolation_by_distance = TRUE
) {
  # Initialize distance matrix
  dist_matrix <- matrix(0, nrow = n_pops, ncol = n_pops)
  
  for (i in 1:n_pops) {
    for (j in i:n_pops) {
      if (i != j) {
        # Calculate geographical Euclidean distance
        geo_dist <- sqrt(sum((pop_positions[i, 1:geo_dims] - pop_positions[j, 1:geo_dims])^2))
        
        # Calculate genetic Euclidean distance if genetic dimensions exist
        gen_dist <- 0
        if (ncol(pop_positions) > geo_dims) {
          gen_dims <- (geo_dims+1):ncol(pop_positions)
          gen_dist <- sqrt(sum((pop_positions[i, gen_dims] - pop_positions[j, gen_dims])^2))
        }
        
        # Combine geographic and genetic distances
        if (use_isolation_by_distance) {
          # With isolation by distance, weight geographic distance more
          raw_dist <- isolation_factor * geo_dist + (1 - isolation_factor) * gen_dist
        } else {
          # Without isolation by distance, use complete distance
          raw_dist <- sqrt(sum((pop_positions[i,] - pop_positions[j,])^2))
        }
        
        # Store the distance
        dist_matrix[i, j] <- raw_dist
        dist_matrix[j, i] <- raw_dist  # Ensure symmetry
      }
    }
  }
  
  return(dist_matrix)
}

#' Check properties of a distance matrix
#' return Prints diagnostics about the matrix
check_distance_matrix <- function(dist_matrix, message = "Distance matrix check", verbose = TRUE) {
  if (!verbose) return()
  
  cat("\n", message, ":\n")
  cat("Dimensions:", dim(dist_matrix), "\n")
  cat("Is symmetric:", all(abs(dist_matrix - t(dist_matrix)) < 1e-10), "\n")
  cat("Range of distances:", range(dist_matrix), "\n")
  cat("Mean distance:", mean(dist_matrix), "\n")
  
  # Check if all values are the same
  if (sd(as.vector(dist_matrix)) < 1e-10) {
    cat("WARNING: All distances are essentially identical!\n")
  }
  
  # Check diagonal
  diag_vals <- diag(dist_matrix)
  cat("Diagonal values - min:", min(diag_vals), "max:", max(diag_vals), "\n")
  
  # Check for skewed distribution
  quants <- quantile(as.vector(dist_matrix))
  cat("Quantiles: 0%:", quants[1], "25%:", quants[2], "50%:", quants[3], 
      "75%:", quants[4], "100%:", quants[5], "\n")
}

#' Transform raw distances to genetic distances with noise
#' return Final genetic distance matrix
transform_to_genetic_distances <- function(
    raw_distances,
    n_pops,
    nonlinear_factor = 0.7,
    noise_level = 0.1,
    use_nonlinear = TRUE,
    use_noise = TRUE,
    verbose = TRUE
) {
  # Check the raw distances
  check_distance_matrix(raw_distances, "Raw distances before transformation", verbose)
  
  # # If raw distances are all very small or very similar, scale them to ensure variation
  # if (max(raw_distances) < 0.1 || sd(as.vector(raw_distances)) < 0.01) {
  #   if (verbose) cat("WARNING: Raw distances have very little variation. Applying scaling.\n")
  #   raw_distances <- raw_distances * (5 / max(raw_distances))
  # }
  
  # Initialize genetic distance matrix
  genetic_dist <- matrix(0, nrow = n_pops, ncol = n_pops)
  
  for (i in 1:n_pops) {
    for (j in i:n_pops) {
      if (i != j) {
        raw_dist <- raw_distances[i, j]
        
        # Apply non-linear transformation if requested
        if (use_nonlinear) {
          # Similar to Fst transformation: Fst = dist/(1 + dist)
          dist_val <- raw_dist^nonlinear_factor / (1 + raw_dist^nonlinear_factor)
        } else {
          # Linear scaling to 0-1 range
          dist_val <- raw_dist / max(raw_distances)
        }
        
        # Add structured noise if requested
        if (use_noise) {
          # More variation for larger distances
          noise_factor <- noise_level * (0.5 + dist_val)
          dist_val <- dist_val * (1 + rnorm(1, 0, noise_factor))
        }
        
        # Ensure distance is non-negative and properly bounded (0-1 range for FST-like measure)
        dist_val <- max(0, min(1, dist_val))
        
        # Store the distance
        genetic_dist[i, j] <- dist_val
        genetic_dist[j, i] <- dist_val  # Ensure symmetry
      }
    }
  }
  
  # Check the final distances
  check_distance_matrix(genetic_dist, "Final genetic distances after transformation", verbose)
  
  return(genetic_dist)
}

# 4. VISUALIZATION FUNCTIONS ####################################################

#' Create a heatmap of genetic distances
#' return heatmap object without plotting it
create_distance_heatmap <- function(dist_matrix, pop_info) {
  # Check if required packages are installed
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    install.packages("ComplexHeatmap")
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    install.packages("circlize")
  }
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    install.packages("RColorBrewer")
  }
  
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  
  # Make sure to assign proper row and column names to the distance matrix
  rownames(dist_matrix) <- pop_info$Population
  colnames(dist_matrix) <- pop_info$Population
  
  # Create annotation dataframe for the heatmap
  annotation_df <- data.frame(
    MajorGroup = pop_info$MajorGroup,
    row.names = pop_info$Population
  )
  
  # Define colors for the major groups
  n_groups <- length(unique(pop_info$MajorGroup))
  group_colors <- brewer.pal(min(n_groups, 9), "Set1")
  if (n_groups > 9) {
    # Add more colors if needed
    group_colors <- c(group_colors, brewer.pal(min(n_groups - 9, 8), "Set2"))
  }
  
  # Create annotation colors
  group_col <- setNames(group_colors[1:n_groups], unique(pop_info$MajorGroup))
  
  # Create row and column annotations
  ha_row <- ComplexHeatmap::rowAnnotation(
    MajorGroup = annotation_df$MajorGroup,
    col = list(MajorGroup = group_col),
    show_legend = TRUE
  )
  
  ha_col <- ComplexHeatmap::HeatmapAnnotation(
    MajorGroup = annotation_df$MajorGroup,
    col = list(MajorGroup = group_col),
    show_legend = FALSE
  )
  
  # Create heatmap without displaying it
  max_dist <- max(dist_matrix)
  heatmap_obj <- ComplexHeatmap::Heatmap(
    dist_matrix,
    name = "GeneticDist",
    col = colorRamp2(c(0, max_dist), c("white", "#1f77b4")),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = (nrow(dist_matrix) <= 50), # Show names only for smaller matrices
    show_column_names = (nrow(dist_matrix) <= 50),
    row_names_side = "left",
    column_names_side = "top",
    row_names_gp = gpar(fontsize = 9),
    column_names_gp = gpar(fontsize = 9),
    show_heatmap_legend = TRUE,
    left_annotation = ha_row,
    top_annotation = ha_col,
    heatmap_legend_param = list(
      title = "Genetic Distance", 
      title_position = "leftcenter-rot"
    )
  )
  
  return(heatmap_obj)
}

#' Create MDS plot of genetic distances
#' return ggplot object with MDS plot without displaying it
create_mds_plot <- function(dist_matrix, pop_info) {
  # Check if required packages are installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
  }
  
  library(ggplot2)
  
  # Make sure to assign proper row names to the distance matrix
  rownames(dist_matrix) <- pop_info$Population
  colnames(dist_matrix) <- pop_info$Population
  
  # For MDS, we need to convert our distance matrix to a dist object
  dist_obj <- as.dist(dist_matrix)
  
  # Calculate MDS coordinates using cmdscale (part of base R)
  mds_coords <- as.data.frame(cmdscale(dist_obj, k = 2))
  colnames(mds_coords) <- c("Dim1", "Dim2")
  
  # Add population information
  mds_coords$Population <- rownames(mds_coords)
  mds_coords$MajorGroup <- pop_info$MajorGroup[match(mds_coords$Population, pop_info$Population)]
  mds_coords$Subgroup <- pop_info$Subgroup[match(mds_coords$Population, pop_info$Population)]
  mds_coords$Admixed <- pop_info$Admixed[match(mds_coords$Population, pop_info$Population)]
  mds_coords$Bottlenecked <- pop_info$Bottlenecked[match(mds_coords$Population, pop_info$Population)]
  
  # Create base MDS plot
  plot <- ggplot(mds_coords, aes(x = Dim1, y = Dim2, color = MajorGroup)) +
    geom_point(aes(shape = ifelse(Bottlenecked == "Yes", "Bottlenecked", "Normal")), 
               size = 3, alpha = 0.8) +
    scale_shape_manual(values = c("Bottlenecked" = 17, "Normal" = 16)) +
    theme_minimal() +
    labs(title = "MDS Plot of Simulated Genetic Distances",
         x = "Dimension 1", y = "Dimension 2",
         shape = "Population Type") +
    theme(legend.position = "right")
  
  # Special highlight for admixed populations
  if (any(pop_info$Admixed == "Yes")) {
    plot <- plot + 
      geom_point(data = subset(mds_coords, Admixed == "Yes"), 
                 size = 5, shape = 1, color = "black")
  }
  
  return(plot)
}

#' Visualize the results of genetic distance simulation
#' return List of plotting functions
visualize_results <- function(sim_results) {
  # Extract data
  dist_matrix <- sim_results$distance_matrix
  pop_info <- sim_results$population_info
  
  # Create visualizations (but don't display them)
  heatmap <- create_distance_heatmap(dist_matrix, pop_info)
  mds_plot <- create_mds_plot(dist_matrix, pop_info)
  
  # Return functions that display plots on new devices
  list(
    heatmap = function() {
      dev.new()
      print(heatmap)
    },
    mds = function() {
      dev.new()
      print(mds_plot)
    }
  )
}

# 5. MAIN SIMULATION FUNCTION ####################################################
# Genetic Distance Separation Parameters
# group_separation and subgroup_separation share the same scale:
#
# Approximate FST equivalent after transformation:
# Values 5-10:  Subtle differentiation    (FST ~0.05-0.15)
# Values 10-20: Moderate differentiation  (FST ~0.15-0.25)
# Values 20-40: Strong differentiation    (FST ~0.25-0.35)
# Values >40:   Extreme differentiation   (FST >0.35)
#
# The parameters differ only in application:
# - group_separation:    Distance between major groups (e.g., continents)
# - subgroup_separation: Distance between subpopulations within major groups
#
# For realistic hierarchical structure, set subgroup_separation to ~1/3 of group_separation
#' Simulate genetic distances using realistic population structure
#' return List containing distance matrix and population information
simulate_genetic_distances <- function(
    n_pops = 50,               # Number of populations
    n_major_groups = 5,        # Number of major population groups (e.g. continental)
    n_subgroups = 12,          # Number of subgroups (e.g. regional)
    geo_dims = 2,              # Geographic dimensions
    genetic_dims = 2,          # Additional genetic drift dimensions
    group_separation = 15,     # Distance between major group centers
    subgroup_separation = 5,   # Distance between subgroup centers
    pop_dispersion = 0.5,      # Dispersion of populations within subgroups
    isolation_factor = 0.8,    # Strength of isolation by distance (0-1)
    admixture_prob = 0.1,      # Proportion of populations with admixture
    bottleneck_prob = 0.05,    # Proportion of populations with bottleneck history
    noise_level = 0.1,         # Level of random noise (measurement error)
    nonlinear_factor = 0.7,    # Non-linearity in distance transformation
    
    # Feature toggle options
    use_subgroups = TRUE,       # Use hierarchical subgroup structure within major groups
    use_genetic_dims = TRUE,    # Use additional genetic dimensions beyond geography
    use_admixture = TRUE,       # Include admixed populations
    use_bottlenecks = TRUE,     # Include bottlenecked populations
    use_isolation_by_distance = TRUE,  # Model isolation by distance
    use_nonlinear = TRUE,       # Use non-linear transformation of distances
    use_noise = TRUE,           # Add stochastic noise to distances
    verbose = TRUE              # Print detailed information
) {
  set.seed(233)  # For reproducibility
  
  # Validate input parameters
  n_pops <- max(10, n_pops)  # Ensure minimum number of populations
  n_major_groups <- min(n_pops/2, max(2, n_major_groups))  # At least 2 groups, at most n_pops/2
  n_subgroups <- min(n_pops, max(n_major_groups, n_subgroups))  # At least as many as major groups
  
  if (verbose) {
    cat("Simulation parameters:\n")
    cat("Populations:", n_pops, "\n")
    cat("Major groups:", n_major_groups, "\n")
    if (use_subgroups) cat("Subgroups:", n_subgroups, "\n")
    cat("Group separation:", group_separation, "\n")
    cat("Subgroup separation:", subgroup_separation, "\n")
    cat("Population dispersion:", pop_dispersion, "\n")
    
    cat("\nEnabled features:\n")
    cat("Subgroups:", use_subgroups, "\n")
    cat("Genetic dimensions:", use_genetic_dims, "\n")
    cat("Admixture:", use_admixture, "\n")
    cat("Bottlenecks:", use_bottlenecks, "\n")
    cat("Isolation by distance:", use_isolation_by_distance, "\n")
    cat("Non-linear transformation:", use_nonlinear, "\n")
    cat("Stochastic noise:", use_noise, "\n")
  }
  
  # Adjust dimensions based on options
  if (!use_genetic_dims) genetic_dims <- 0
  n_dims <- geo_dims + genetic_dims
  
  # Step 1: Create population metadata with group assignments
  pop_info <- create_population_metadata(n_pops, n_major_groups, n_subgroups, use_subgroups)
  
  # Step 2: Generate group and subgroup centroids
  group_centroids <- generate_group_centroids(n_major_groups, n_dims, geo_dims, 
                                              group_separation, use_genetic_dims)
  
  subgroup_centroids <- generate_subgroup_centroids(n_subgroups, n_major_groups, 
                                                    group_centroids, n_dims, 
                                                    subgroup_separation, use_subgroups)
  
  # Step 3: Position populations in multidimensional space
  pop_positions <- position_populations(n_pops, pop_info, subgroup_centroids, 
                                        n_dims, pop_dispersion)
  
  # Step 4: Create populations with special demographic histories
  # First admixed populations (if enabled)
  admix_result <- create_admixed_populations(pop_positions, pop_info, n_pops, 
                                             subgroup_centroids, n_dims, 
                                             admixture_prob, pop_dispersion, use_admixture)
  pop_positions <- admix_result$positions
  pop_info <- admix_result$metadata
  
  # Then bottlenecked populations (if enabled)
  bottleneck_result <- create_bottlenecked_populations(pop_positions, pop_info, n_pops, 
                                                       subgroup_centroids, geo_dims, n_dims, 
                                                       bottleneck_prob, pop_dispersion, use_bottlenecks)
  pop_positions <- bottleneck_result$positions
  pop_info <- bottleneck_result$metadata
  
  # Step 5: Calculate distances
  raw_distances <- calculate_raw_distances(pop_positions, n_pops, geo_dims, 
                                           isolation_factor, use_isolation_by_distance)
  
  genetic_distances <- transform_to_genetic_distances(raw_distances, n_pops, 
                                                      nonlinear_factor, noise_level,
                                                      use_nonlinear, use_noise, verbose)
  
  genetic_distances <- (genetic_distances + t(genetic_distances)) / 2
  rownames(genetic_distances) <- pop_info$Population
  colnames(genetic_distances) <- pop_info$Population
  
  # Return results
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
      use_noise = use_noise
    )
  )
}

# 6.1 USER INTERFACE FUNCTION ####################################################
#' Run a complete genetic distance simulation with visualization
#' @param model Character string specifying the simulation model:
#'   "mixed" - Both geographic and genetic factors (default)
#'   "geography" - Geography-only model
#'   "genetics" - Genetics-only model
#'   "custom" - Custom parameter configuration
#' return List containing simulation results and plots
run_genetic_simulation <- function(
    # Basic parameters
  n_pops = 30,                 # Number of populations
  n_major_groups = 4,          # Number of major groups
  n_subgroups = 8,             # Number of subgroups
  
  # Model selection
  model = "mixed",             # "mixed", "geography", "genetics", or "custom"
  
  # Geography parameters
  geo_dims = NULL,             # Geographic dimensions (auto-set based on model)
  isolation_factor = NULL,     # Geography-genetics balance (auto-set based on model)
  
  # Genetic parameters
  genetic_dims = NULL,         # Genetic dimensions (auto-set based on model)
  
  # Structural parameters
  group_separation = 15,       # Spatial separation between groups
  subgroup_separation = NULL,  # Distance between subgroups (default: group_separation/3)
  pop_dispersion = 0.5,        # Within-population dispersion
  
  # Demographic history
  admixture_prob = 0.15,       # Proportion of admixed populations
  bottleneck_prob = 0.10,      # Proportion of bottlenecked populations
  
  # Feature toggles (automatically adjusted based on model)
  use_subgroups = TRUE,        # Use hierarchical subgroup structure
  use_genetic_dims = NULL,     # Use additional genetic dimensions
  use_admixture = TRUE,        # Include admixed populations
  use_bottlenecks = TRUE,      # Include bottlenecked populations
  use_isolation_by_distance = NULL,  # Model isolation by distance
  use_nonlinear = TRUE,        # Use non-linear transformation
  use_noise = TRUE,            # Add stochastic noise
  
  # Output options
  output_file = NULL,          # Optional file to save results
  verbose = TRUE               # Print diagnostic information
) {
  # Implement model-specific parameter configurations
  if (model == "geography") {
    # Geography-only model
    geo_dims <- ifelse(is.null(geo_dims), 2, geo_dims)
    genetic_dims <- 0
    isolation_factor <- 1.0
    use_isolation_by_distance <- TRUE
    use_genetic_dims <- FALSE
    
  } else if (model == "genetics") {
    # Genetics-only model
    geo_dims <- 0
    genetic_dims <- ifelse(is.null(genetic_dims), 2, genetic_dims)
    isolation_factor <- 0.0
    use_isolation_by_distance <- FALSE
    use_genetic_dims <- TRUE
    
  } else if (model == "mixed") {
    # Mixed model (default)
    geo_dims <- ifelse(is.null(geo_dims), 2, geo_dims)
    genetic_dims <- ifelse(is.null(genetic_dims), 2, genetic_dims)
    isolation_factor <- ifelse(is.null(isolation_factor), 0.8, isolation_factor)
    use_isolation_by_distance <- TRUE
    use_genetic_dims <- TRUE
    
  } else if (model == "custom") {
    # Custom model - use provided parameters or defaults
    geo_dims <- ifelse(is.null(geo_dims), 2, geo_dims)
    genetic_dims <- ifelse(is.null(genetic_dims), 2, genetic_dims)
    isolation_factor <- ifelse(is.null(isolation_factor), 0.8, isolation_factor)
    use_isolation_by_distance <- ifelse(is.null(use_isolation_by_distance), 
                                        isolation_factor > 0, use_isolation_by_distance)
    use_genetic_dims <- ifelse(is.null(use_genetic_dims), 
                               genetic_dims > 0, use_genetic_dims)
  } else {
    stop("Invalid model specified. Choose from: 'mixed', 'geography', 'genetics', or 'custom'")
  }
  
  # Ensure parameters are internally consistent
  if (use_isolation_by_distance && is.null(isolation_factor)) {
    isolation_factor <- 0.8  # Default if isolation by distance is enabled
  }
  
  if (isolation_factor == 1.0) {
    use_isolation_by_distance <- TRUE
  }
  
  if (geo_dims == 0) {
    use_isolation_by_distance <- FALSE
    isolation_factor <- 0
  }
  
  if (genetic_dims == 0) {
    use_genetic_dims <- FALSE
  }
  
  if (use_genetic_dims && genetic_dims == 0) {
    genetic_dims <- 2  # Default if genetic dimensions are enabled
  }
  
  # Set default subgroup separation if not provided
  if (is.null(subgroup_separation)) {
    subgroup_separation <- group_separation / 3
  }
  
  # Print model information
  if (verbose) {
    cat("Running simulation with", model, "model\n")
    cat("Geographic dimensions:", geo_dims, "\n")
    cat("Genetic dimensions:", genetic_dims, "\n")
    if (use_isolation_by_distance) {
      cat("Isolation factor:", isolation_factor, 
          "(higher values give more weight to geography)\n")
    }
  }
  
  # Run simulation with the configured parameters
  if (verbose) cat("Simulating genetic distances...\n")
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
    verbose = verbose
  )
  
  # Generate visualizations
  if (verbose) cat("Creating visualizations...\n")
  plots <- visualize_results(sim_results)
  
  # Print simulation summary
  if (verbose) {
    cat("\nSimulation Summary:\n")
    cat("Model type:", model, "\n")
    cat("Number of populations:", n_pops, "\n")
    cat("Number of major groups:", n_major_groups, "\n")
    if (use_subgroups) cat("Number of subgroups:", n_subgroups, "\n")
    cat("Geographic dimensions:", geo_dims, "\n")
    cat("Genetic dimensions:", genetic_dims, "\n")
    if (use_isolation_by_distance) cat("Isolation factor:", isolation_factor, "\n")
    cat("Group separation:", group_separation, "\n")
    cat("Subgroup separation:", subgroup_separation, "\n")
    if (use_admixture) cat("Admixed populations:", sum(sim_results$population_info$Admixed == "Yes"), "\n")
    if (use_bottlenecks) cat("Bottlenecked populations:", sum(sim_results$population_info$Bottlenecked == "Yes"), "\n")
  }
  
  # Save results to file if requested
  if (!is.null(output_file)) {
    # Add .RData extension if not present
    if (!grepl("\\.RData$", output_file)) {
      output_file <- paste0(output_file, ".RData")
    }
    
    if (verbose) cat("Saving results to", output_file, "\n")
    save(sim_results, plots, file = output_file)
    
    # Also export CSV files for easier access
    csv_output <- gsub("\\.RData$", ".csv", output_file)
    write.csv(sim_results$distance_matrix, csv_output)
  }
  
  # Return results and plots
  if (verbose) cat("\nReturning results and plots. Use print() to display them.\n")
  return(list(
    results = sim_results,
    plots = plots
  ))
}

# 6.2 PREDEFINED SCENARIOS ####################################################

#' Run simulation with predefined biological scenarios
#' return Simulation results for the selected scenario
run_genetic_scenario <- function(
    scenario = "default",      # Scenario name
    n_pops = 30,               # Number of populations to simulate
    output_file = NULL,        # Optional file to save results
    verbose = TRUE             # Print diagnostic information
) {
  # Define parameter sets for different biological scenarios
  scenarios <- list(
    # Basic simulation with all features enabled
    "default" = list(
      model = "mixed",
      n_major_groups = 4,
      n_subgroups = 8,
      admixture_prob = 0.15,
      bottleneck_prob = 0.1,
      group_separation = 15,
      subgroup_separation = 5,
      pop_dispersion = 0.5,
      use_subgroups = TRUE,
      use_admixture = TRUE,
      use_bottlenecks = TRUE,
      use_nonlinear = TRUE,
      use_noise = TRUE
    ),
    
    # Island model: Well-separated populations with minimal gene flow
    "island" = list(
      model = "mixed",
      n_major_groups = 6,
      n_subgroups = 6,  # One subgroup per major group
      admixture_prob = 0.05,  # Very low admixture
      bottleneck_prob = 0.2,  # High bottleneck probability (founder effects)
      group_separation = 25,  # High separation between groups
      subgroup_separation = 0,  # No subgroup structure
      pop_dispersion = 0.3,   # Tight clusters (low within-island diversity)
      use_subgroups = FALSE,  # No subgroups (each island is distinct)
      use_admixture = TRUE,
      use_bottlenecks = TRUE,
      use_nonlinear = TRUE,
      use_noise = TRUE
    ),
    
    # Stepping stone: Continuous distribution with isolation by distance
    "stepping_stone" = list(
      model = "geography",  # Geography dominates
      n_major_groups = 3,
      n_subgroups = 12,  # Many subgroups for finer geographic structure
      admixture_prob = 0.1,
      bottleneck_prob = 0.05,
      group_separation = 10,  # Lower separation
      subgroup_separation = 3,  # Closer subgroups
      pop_dispersion = 0.7,   # More dispersed within groups
      use_subgroups = TRUE,
      use_admixture = TRUE,
      use_bottlenecks = FALSE,
      use_nonlinear = TRUE,
      use_noise = TRUE
    ),
    
    # Recent admixture: Populations with extensive mixing
    "admixture" = list(
      model = "mixed",
      n_major_groups = 3,
      n_subgroups = 6,
      admixture_prob = 0.4,  # High admixture
      bottleneck_prob = 0,   # No bottlenecks
      group_separation = 12,
      subgroup_separation = 4,
      pop_dispersion = 0.6,  # Higher within-group diversity
      use_subgroups = TRUE,
      use_admixture = TRUE,
      use_bottlenecks = FALSE,
      use_nonlinear = TRUE,
      use_noise = TRUE
    ),
    
    # Ancient divergence: Old, well-differentiated populations
    "ancient_divergence" = list(
      model = "genetics",  # Genetics dominates (not geography)
      n_major_groups = 4,
      n_subgroups = 10,
      admixture_prob = 0.05,  # Low admixture
      bottleneck_prob = 0.1,
      group_separation = 30,  # Very high separation
      subgroup_separation = 8,
      pop_dispersion = 0.4,   # Moderate within-group diversity
      use_subgroups = TRUE,
      use_admixture = TRUE,
      use_bottlenecks = TRUE,
      use_nonlinear = TRUE,
      use_noise = FALSE  # Less noise (well-established differences)
    ),
    
    # Simplified model with just major groups
    "simple" = list(
      model = "custom",
      geo_dims = 2,
      genetic_dims = 0,
      n_major_groups = 3,
      n_subgroups = 0,
      admixture_prob = 0,
      bottleneck_prob = 0,
      group_separation = 15,
      subgroup_separation = 0,
      pop_dispersion = 0.5,
      use_subgroups = FALSE,
      use_admixture = FALSE,
      use_bottlenecks = FALSE,
      use_nonlinear = FALSE,
      use_noise = FALSE
    )
  )
  
  # Check if the requested scenario exists
  if (!scenario %in% names(scenarios)) {
    stop("Unknown scenario: ", scenario, 
         ". Available scenarios: ", paste(names(scenarios), collapse = ", "))
  }
  
  # Get parameters for the selected scenario
  params <- scenarios[[scenario]]
  
  # Add scenario name to output file
  if (!is.null(output_file)) {
    # Just append the scenario name to the base filename
    output_file <- paste0(output_file, "_", scenario)
    # The .RData extension will be added by run_genetic_simulation()
  }
  
  # Print scenario information
  if (verbose) {
    cat("Running '", scenario, "' scenario\n", sep = "")
    cat("This models: ", switch(scenario,
                                "default" = "A standard population genetic structure with all features",
                                "island" = "Island model with well-separated populations and founder effects",
                                "stepping_stone" = "Continuous populations with isolation by distance",
                                "admixture" = "Populations with extensive recent admixture",
                                "ancient_divergence" = "Deeply diverged populations with little gene flow",
                                "simple" = "Simplified model with just major groups"), "\n", sep = "")
  }
  
  # Run the simulation with the selected parameters and verbose setting
  do.call(run_genetic_simulation, c(list(n_pops = n_pops, 
                                         output_file = output_file,
                                         verbose = verbose), 
                                    params))
}

# USAGE ####################################################

# Run a simulation with a specific scenario
if (F) {
  result <- run_genetic_scenario(
    scenario = "island",         # Choose a predefined scenario
    n_pops = 30,                 # Number of populations to simulate
    output_file = "./genetic_simulation_results"  # Optional: save results to file
  )
}

# Genetics-based simulation (no geographic factors)
if (T) {
  result <- run_genetic_simulation(
    model = "genetics",
    n_pops = 30,                  # Number of populations
    n_major_groups = 3,           # Number of major groups
    n_subgroups = 8,              # Number of subgroups
    admixture_prob = 0.15,        # Proportion of admixed populations
    bottleneck_prob = 0.1,        # Proportion of bottlenecked populations
    group_separation = 3,         # Distance between groups
    subgroup_separation = 1,      # Distance between subgroups
    pop_dispersion = 0.5,         # Within-population dispersion
    
    output_file = "./GDM_simulated",  # Optional: save results to file
  )
}

# run with custom settings
if (F) {
  result <- run_genetic_simulation(
    model = "custom",             # Using custom configuration
    n_pops = 30,                  # Number of populations
    n_major_groups = 3,           # Number of major groups
    n_subgroups = 8,              # Number of subgroups
    admixture_prob = 0.15,        # Proportion of admixed populations
    bottleneck_prob = 0.1,        # Proportion of bottlenecked populations
    group_separation = 15,        # Distance between groups
    subgroup_separation = 5,      # Distance between subgroups
    pop_dispersion = 0.5,         # Within-population dispersion
    
    # Parameter configuration
    geo_dims = 2,                 # Using 2D geography
    genetic_dims = 0,             # No genetic dimensions
    isolation_factor = 0,         # No isolation by distance effect
    
    # Feature toggles
    use_subgroups = TRUE,         # Use hierarchical subgroup structure
    use_genetic_dims = FALSE,     # No additional genetic dimensions
    use_admixture = FALSE,        # No admixed populations
    use_bottlenecks = FALSE,      # No bottlenecked populations
    use_isolation_by_distance = FALSE,  # No isolation by distance
    use_nonlinear = FALSE,        # Linear transformation
    use_noise = TRUE              # Add stochastic noise
  )
}

# display the plots
# result$plots$heatmap()
# result$plots$mds()

# Access the genetic distance matrix
# dist_matrix <- result$results$distance_matrix
# View(dist_matrix)

# Access population metadata
# pop_info <- result$results$population_info
# View(pop_info)
