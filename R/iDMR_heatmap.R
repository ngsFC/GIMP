#' Generate Heatmap of Imprinted DMRs Methylation
#'
#' This function generates a heatmap for visualizing methylation data of Imprinted Differentially Methylated Regions (DMRs).
#'
#' @param df_ICR A data frame or matrix containing methylation beta values for Imprinted DMRs.
#' @param group_vector A vector indicating the group labels (e.g., "Control" and "Case") for each sample in `df_ICR`.
#' Each element in `group_vector` should correspond to a sample in `df_ICR`.
#' @param control_label A character string specifying the label for the control group in `group_vector`. Default is `"Control"`.
#' @param case_label A character string specifying the label for the case group in `group_vector`. Default is `"Case"`.
#' @param bedmeth A character string specifying the BED data version for DMR coordinates. Options are `"v1"`, `"v2"`, or `"450k"`. Default is `"v1"`.
#' @param order_by A character string specifying the ordering rows in the heatmap. Options are `"cord"` for coordinates or `"meth"` for methylation values. Default is `"cord"`.
#' @param annotation_col A named list of colors for each unique value in `group_vector`. If `NULL`, default colors are assigned using the "viridis" palette. Default is `NULL`.
#' @param plot_type A character string specifying the type of heatmap to generate. Options are `"beta"` for beta values, `"delta"` for values normalized against controls, and `"defect"` for defect matrix based on standard deviations. Default is `"beta"`.
#' @param sd_threshold A numeric value specifying the standard deviation threshold for detecting defects in the defect matrix. Only used if `plot_type` is `"defect"`. Default is `3`.
#' @return A heatmap plot visualizing methylation of Imprinted DMRs.
#' @import pheatmap
#' @import viridisLite
#' @import ggplotify
#' @examples
#' # Example group_vector with "Case" and "Control" labels for each sample
#' group_vector <- c(rep("Case", 10), rep("Control", 10))
#' DMR_heatmap(df_ICR = my_ICR_data, group_vector = group_vector, annotation_col = list(Sample = c("darkgreen", "darkred")))
#' @export

iDMR_heatmap <- function(df_ICR, group_vector, control_label = "Control", case_label = "Case", bedmeth = "v1", order_by = "cord", annotation_col = NULL, plot_type = "beta", sd_threshold = 3) {
  
  # Load BED data based on bedmeth version
  if (bedmeth == "v1" || bedmeth == "450k") {
    data(DMRs.hg19)
    ICR_cord <- DMRs.hg19
    odr <- ICR_cord$ICR
  } else if (bedmeth == "v2") {
    data(DMRs.hg38)
    ICR_cord <- DMRs.hg38
    odr <- ICR_cord$ICR
  } else {
    stop("Invalid bedmeth version. Choose from 'v1', 'v2', or '450k'.")
  }
  
  # Ensure group_vector length matches number of samples in input data
  if (length(group_vector) != ncol(df_ICR)) {
    stop("Length of 'group_vector' must match the number of samples (columns) in 'df_ICR'.")
  }
  
  group_vector <- factor(group_vector, levels = c(control_label, case_label))
  annotation_name <- if (!is.null(annotation_col)) names(annotation_col)[1] else "Sample"
  
  # Preparing annotation data
  mat_col <- data.frame(group_vector)
  colnames(mat_col) <- annotation_name
  rownames(mat_col) <- colnames(df_ICR)
  
  unique_groups <- levels(group_vector)
  if (is.null(annotation_col)) {
    default_colors <- viridisLite::viridis(length(unique_groups))
    annotation_col <- list(setNames(default_colors, unique_groups))
    names(annotation_col) <- annotation_name
  } else {
    if (!is.list(annotation_col) || length(annotation_col[[1]]) != length(unique_groups)) {
      stop("The 'annotation_col' list must have the same number of colors as the unique values in 'group_vector'.")
    }
    names(annotation_col[[1]]) <- unique_groups
  }
  
  # Define heatmap palette and breaks
  paletteLength <- 100
  
  if (plot_type == "beta") {
    # Define palette and breaks for beta values
    colorPalette <- colorRampPalette(c("#785EF0", "white", "#9a031e"))(paletteLength)
    breaks <- c(seq(0, 0.5, length.out = ceiling(paletteLength / 2) + 1), 
                seq(0.500001, 1, length.out = floor(paletteLength / 2)))
    heatmap_data <- df_ICR[odr, ]
    main_title <- "Methylation of Imprinted DMRs (Beta)"
    
  } else if (plot_type == "delta") {
    # Calculate DeltaBeta matrix
    control_indices <- which(group_vector == control_label)
    control_means <- rowMeans(df_ICR[, control_indices, drop = FALSE])
    heatmap_data <- sweep(df_ICR, 1, control_means)
    colorPalette <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
    breaks <- seq(-0.5, 0.5, length.out = paletteLength + 1)
    main_title <- "Delta Beta of Imprinted DMRs"
    
  } else if (plot_type == "defect") {
    # Calculate binary defect matrix based on sd_threshold
    df_mvalues <- log2(df_ICR / (1 - df_ICR))
    control_indices <- which(group_vector == control_label)
    control_mvalues <- df_mvalues[, control_indices, drop = FALSE]
    control_means_m <- rowMeans(control_mvalues)
    control_sds_m <- apply(control_mvalues, 1, sd)
    upper_threshold <- control_means_m + sd_threshold * control_sds_m
    lower_threshold <- control_means_m - sd_threshold * control_sds_m
    heatmap_data <- (df_mvalues > upper_threshold | df_mvalues < lower_threshold) * 1
    colorPalette <- c("white", "black")
    breaks <- NULL  # No breaks needed for binary data
    main_title <- "Defect Matrix of Imprinted DMRs"
    
  } else {
    stop("Invalid plot_type. Choose from 'beta', 'delta', or 'defect'.")
  }
  
  # Ordering rows
  row_clust <- if (order_by == "cord") FALSE else if (order_by == "meth") TRUE else {
    stop("Please select a valid 'order_by' parameter: 'cord' to order by coordinates or 'meth' to order by methylation values.")
  }
  
  # Generate the heatmap
  pheatmap(
    mat = heatmap_data,
    annotation_col = mat_col,
    annotation_colors = annotation_col,
    color = colorPalette,
    breaks = breaks,
    border_color = "grey",
    main = main_title,
    annotation_legend = TRUE,
    annotation_names_col = FALSE,
    annotation_names_row = FALSE,
    drop_levels = FALSE,
    fontsize = 8,
    cluster_rows = row_clust,
    cluster_cols = TRUE,
    clustering_distance_rows = if (plot_type == "defect") "binary" else "euclidean",
    clustering_distance_cols = "euclidean", 
    clustering_method = "ward.D2"
  ) %>% as.ggplot()
}
