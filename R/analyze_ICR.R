#' Create the ICR CPG Matrix
#'
#' @examples
#'  Assuming df_ICR is your matrix and you have a group_vector like c("Control", "Case", "Control", "Case", ...)
#'  result <- analyze_ICR(df_ICR, group_vector)
#'  View the heatmaps and matrices:
#'  result$heatmap_ICR  # Heatmap of df.ICR
#'  result$df_ICR_delta  # DeltaBeta matrix
#'  result$df_ICR_defect  # Binary defect matrix
#' @export


# Function to create heatmaps and delta matrices
analyze_ICR <- function(df_ICR, group_vector, control_label = "Control", case_label = "Case") {
  
  # Ensure the group_vector is a factor
  group_vector <- factor(group_vector, levels = c(control_label, case_label))
  
  # Step 1: Plot a heatmap of df.ICR with Case/Control annotations
  plot_heatmap <- function(df, group_vector, title = "Heatmap of df.ICR") {
    df_melt <- melt(as.matrix(df))
    colnames(df_melt) <- c("ICR", "Sample", "BetaValue")
    
    # Add the group annotations (Case/Control)
    df_melt$Group <- group_vector[df_melt$Sample]
    
    # Plot the heatmap
    ggplot(df_melt, aes(x = Sample, y = ICR, fill = BetaValue)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.5) +
      labs(title = title, x = "Sample", y = "ICR", fill = "Beta Value") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
  
  # Step 2: Calculate DeltaBeta (df.ICR.delta)
  calculate_delta_beta <- function(df, group_vector, control_label) {
    control_indices <- which(group_vector == control_label)
    
    # Calculate the mean for controls for each ICR (row-wise)
    control_means <- rowMeans(df[, control_indices, drop = FALSE])
    
    # Subtract control means from the entire matrix
    df_delta <- sweep(df, 1, control_means)
    
    return(df_delta)
  }
  
  # Step 3: Create a binary defect matrix (df.ICR.defect)
  calculate_defect_matrix <- function(df, group_vector, control_label) {
    control_indices <- which(group_vector == control_label)
    
    # Calculate the mean and 3*SD for controls for each ICR (row-wise)
    control_means <- rowMeans(df[, control_indices, drop = FALSE])
    control_sds <- apply(df[, control_indices, drop = FALSE], 1, sd)
    threshold <- control_means + 3 * control_sds
    
    # Create a binary matrix where BetaValue > threshold is 1, otherwise 0
    df_defect <- sweep(df, 1, threshold, FUN = ">") * 1  # Binary matrix (1 if exceeds, else 0)
    
    return(df_defect)
  }
  
  # Step 1: Heatmap of df.ICR with Case/Control annotations
  print("Plotting heatmap of df.ICR")
  heatmap_1 <- plot_heatmap(df_ICR, group_vector, title = "Heatmap of df.ICR")
  
  # Step 2: Calculate DeltaBeta matrix and plot heatmap
  df_ICR_delta <- calculate_delta_beta(df_ICR, group_vector, control_label)
  print("Plotting heatmap of DeltaBeta (df_ICR_delta)")
  heatmap_2 <- plot_heatmap(df_ICR_delta, group_vector, title = "Heatmap of DeltaBeta (df_ICR_delta)")
  
  # Step 3: Calculate binary defect matrix and plot heatmap
  df_ICR_defect <- calculate_defect_matrix(df_ICR, group_vector, control_label)
  print("Plotting heatmap of binary defect matrix (df_ICR_defect)")
  heatmap_3 <- plot_heatmap(df_ICR_defect, group_vector, title = "Heatmap of Defect Matrix (df_ICR_defect)")
  
  # Return the three matrices and their heatmaps
  return(list(
    df_ICR_delta = df_ICR_delta,
    df_ICR_defect = df_ICR_defect,
    heatmap_ICR = heatmap_1,
    heatmap_delta = heatmap_2,
    heatmap_defect = heatmap_3
  ))
}

