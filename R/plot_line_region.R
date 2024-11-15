plot_line_region <- function(significantDMPs, ICRcpg, ICR, sampleInfo) {

regionDMPs <- significantDMPs[significantDMPs$ICR == ICR, ]
  if (nrow(regionDMPs) == 0) {
    stop("No significant DMPs found for the specified ICR.")
  }
  
  regionStart <- min(regionDMPs$start)
  regionEnd <- max(regionDMPs$end)
  
  regionCpGs <- ICRcpg[ICRcpg$cstart >= regionStart & ICRcpg$cstart <= regionEnd, ]
  if (nrow(regionCpGs) == 0) {
    stop("No CpGs found in the specified ICR region.")
  }
  
  methylationData <- regionCpGs[, 1:(ncol(regionCpGs) - 4)]
  annotationData <- regionCpGs[, (ncol(regionCpGs) - 3):ncol(regionCpGs)]
  
  methylationData$cstart <- annotationData$cstart
  methylationDataLong <- reshape2::melt(methylationData, id.vars = "cstart", variable.name = "Sample", value.name = "Methylation")
  methylationDataLong$Type <- rep(sampleInfo, each = nrow(methylationData))
  methylationDataLong$IsDMP <- methylationDataLong$cstart %in% regionDMPs$cstart
  
  plot <- ggplot(methylationDataLong, aes(x = cstart, y = Methylation, group = Sample, color = Type)) +
    geom_line(aes(linetype = Type), alpha = 0.5) +  # Adjusted alpha for lines to enhance comparability
    geom_point(data = methylationDataLong[!methylationDataLong$IsDMP, ], size = 0.5, color = "grey") +  
    geom_point(data = methylationDataLong[methylationDataLong$IsDMP, ], size = 0.5, color = "black") + 
    geom_vline(data = data.frame(x = regionDMPs$cstart), aes(xintercept = x), color = "lightgrey", linetype = "dashed") +
    scale_color_manual(values = c("Control" = "grey", "Case" = "red")) +
    scale_linetype_manual(values = c("Control" = "solid", "Case" = "solid")) +
    labs(
      title = ICR,
      x = "CpG Coordinates",
      y = "Methylation Value"
    ) +
    theme_minimal()
  
  return(plot)
}

