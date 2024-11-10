#' Generate Heatmap of Imprinted DMRs Methylation
#'
#' This function generates a heatmap for visualizing methylation data of Imprinted Differentially Methylated Regions (DMRs).
#' It allows custom color schemes for group annotations, with default colors from the "viridis" palette.
#'
#' @param df_ICR A data frame or matrix containing methylation beta values for Imprinted DMRs.
#' @param group_vector A vector indicating the group labels (e.g., "Control" and "Case") for each sample in `df_ICR`.
#' Each element in `group_vector` should correspond to a sample in `df_ICR`.
#' @param control_label A character string specifying the label for the control group in `group_vector`. Default is `"Control"`.
#' @param case_label A character string specifying the label for the case group in `group_vector`. Default is `"Case"`.
#' @param bedmeth A character string specifying the BED data version for DMR coordinates. Options are `"v1"`, `"v2"`, or `"450k"`. Default is `"v1"`.
#' @param cluster_by A character string specifying the clustering method for rows in the heatmap. Options are `"cord"` or `"meth"`. Default is `"cord"`.
#' @param annotation_col A named list of colors for each unique value in `group_vector`. If `NULL`, default colors are assigned using the "viridis" palette. Default is `NULL`.
#' @return A heatmap plot visualizing methylation of Imprinted DMRs.
#' @examples
#' # Example group_vector with "Case" and "Control" labels for each sample
#' group_vector <- c(rep("Case", 10), rep("Control", 10))
#' DMR_heatmap(df_ICR = my_ICR_data, group_vector = group_vector, annotation_col = list(Sample = c("darkgreen", "darkred")))
#' @export


iDMR_heatmap <- function(df_ICR, group_vector, control_label = "Control", case_label = "Case", 
                         annotation_col = NULL, cluster_rows = TRUE, cluster_cols = TRUE) {

library(ggplot2)
library(reshape2)
library(dplyr)
library(stats)
  
  # Ensure group_vector length matches number of samples (columns) in df_ICR
  if (length(group_vector) != ncol(df_ICR)) {
    stop("Length of 'group_vector' must match the number of samples (columns) in 'df_ICR'.")
  }
  
  # Create annotation data frame
  group_vector <- factor(group_vector, levels = c(control_label, case_label))
  annotation_data <- data.frame(Sample = colnames(df_ICR), Group = group_vector)
  rownames(annotation_data) <- annotation_data$Sample
  
  # Define color palettes
  beta_palette <- colorRampPalette(c("#785EF0", "white", "#9a031e"))(100)
  delta_palette <- colorRampPalette(c("blue", "white", "red"))(100)
  binary_palette <- c("white", "black")
  
  create_heatmap <- function(matrix_data, group_vector, title, color_palette, scale_limits = NULL) {
    # Perform clustering if required
    if (cluster_rows) {
      row_dend <- hclust(dist(matrix_data), method = "ward.D2")
      matrix_data <- matrix_data[order.dendrogram(as.dendrogram(row_dend)), , drop = FALSE]
    }
    
    if (cluster_cols) {
      col_dend <- hclust(dist(t(matrix_data)), method = "ward.D2")
      matrix_data <- matrix_data[, order.dendrogram(as.dendrogram(col_dend)), drop = FALSE]
    }
    
    df_melt <- melt(matrix_data)
    colnames(df_melt) <- c("ICR", "Sample", "Value")
    df_melt <- df_melt %>%
      mutate(Group = group_vector[Sample])
    
    ggplot(df_melt, aes(x = Sample, y = ICR, fill = Value)) +
      geom_tile() +
      scale_fill_gradientn(colors = color_palette, limits = scale_limits) +
      labs(title = title, x = "Sample", y = "ICR", fill = "Value") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
  
  # Create Beta Heatmap
  heatmap_beta <- create_heatmap(
    matrix_data = df_ICR,
    group_vector = annotation_data$Group,
    title = "Methylation of Imprinted DMRs (Beta)",
    color_palette = beta_palette,
    scale_limits = c(0, 1)
  )
  
  # Calculate DeltaBeta matrix
  control_indices <- which(group_vector == control_label)
  control_means <- rowMeans(df_ICR[, control_indices, drop = FALSE])
  df_ICR_delta <- sweep(df_ICR, 1, control_means)
  
  # Create Delta Beta Heatmap
  heatmap_delta <- create_heatmap(
    matrix_data = df_ICR_delta,
    group_vector = annotation_data$Group,
    title = "Delta Beta of Imprinted DMRs",
    color_palette = delta_palette,
    scale_limits = c(-0.5, 0.5)
  )
  
  # Calculate binary defect matrix
  df_mvalues <- log2(df_ICR / (1 - df_ICR))
  control_mvalues <- df_mvalues[, control_indices, drop = FALSE]
  control_means_m <- rowMeans(control_mvalues)
  control_sds_m <- apply(control_mvalues, 1, sd)
  upper_threshold <- control_means_m + 3 * control_sds_m
  lower_threshold <- control_means_m - 3 * control_sds_m
  df_ICR_defect <- (df_mvalues > upper_threshold | df_mvalues < lower_threshold) * 1
  
  # Create Defect Matrix Heatmap
  heatmap_defect <- create_heatmap(
    matrix_data = df_ICR_defect,
    group_vector = annotation_data$Group,
    title = "Defect Matrix of Imprinted DMRs",
    color_palette = binary_palette,
    scale_limits = c(0, 1)
  )
  
  # Return list of ggplot heatmaps
  return(list(
    heatmap_beta = heatmap_beta,
    heatmap_delta = heatmap_delta,
    heatmap_defect = heatmap_defect
  ))
}

