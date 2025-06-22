#' Plot Line Plot for ICR Methylation
#'
#' This function generates a line plot to visualize methylation values across a specified ICR.
#' Users can choose between a static `ggplot2` plot or an interactive `plotly` plot.
#'
#' @param significantDMPs A data frame containing information about significant DMPs. Must include columns `ICR`, `start`, and `end`.
#' @param ICRcpg A data frame or matrix containing CpG methylation data. Includes CpG coordinates (`cstart`) and methylation values.
#' @param ICR A character string specifying the name of the ICR region to be plotted.
#' @param sampleInfo A character vector providing group labels (e.g., `"Control"` or `"Case"`) for each sample in the methylation data.
#' @param interactive A logical value indicating whether to return an interactive `plotly` plot (`TRUE`) or a static `ggplot2` plot (`FALSE`). Default is `TRUE`.
#' @return A plot representing the line plot of methylation values across the specified ICR region, highlighting significant DMPs. 
#' The plot is either a `ggplot2` object or a `plotly` object, depending on the value of `interactive`.
#' @examples
#' # Example data for significantDMPs
#' plot <- plot_line_ICR(significantDMPs, ICRcpg, ICR = "KCNQ1OT1:TSS-DMR", sampleInfo = sampleInfo, interactive = T)
#' print(plot)
#' @export

plot_line_ICR <- function(significantDMPs, ICRcpg, ICR, sampleInfo, interactive = TRUE) {
  library(ggplot2)
  library(reshape2)
  if (interactive) library(plotly)
  
  # Filter DMPs for the specified ICR
  regionDMPs <- significantDMPs[significantDMPs$ICR == ICR, ]
  if (nrow(regionDMPs) == 0) {
    stop("No significant DMPs found for the specified ICR.")
  }
  
  regionStart <- min(regionDMPs$start)
  regionEnd <- max(regionDMPs$end)
  
  # Filter CpGs for the specified ICR region
  regionCpGs <- ICRcpg[ICRcpg$cstart >= regionStart & ICRcpg$cstart <= regionEnd, ]
  if (nrow(regionCpGs) == 0) {
    stop("No CpGs found in the specified ICR region.")
  }
  
  # Reshape data for plotting
  methylationData <- regionCpGs[, 1:(ncol(regionCpGs) - 4)]
  annotationData <- regionCpGs[, (ncol(regionCpGs) - 3):ncol(regionCpGs)]
  
  methylationData$cstart <- annotationData$cstart
  methylationDataLong <- reshape2::melt(methylationData, id.vars = "cstart", variable.name = "Sample", value.name = "Methylation")
  methylationDataLong$Type <- rep(sampleInfo, each = nrow(methylationData))
  methylationDataLong$IsDMP <- methylationDataLong$cstart %in% regionDMPs$cstart
  
  plot <- ggplot(methylationDataLong, aes(x = cstart, y = Methylation, group = Sample, color = Type)) +
    geom_line(aes(linetype = Type), alpha = 0.5) +
    geom_point(data = methylationDataLong[!methylationDataLong$IsDMP, ], size = 0.5, color = "grey") +
    geom_point(data = methylationDataLong[methylationDataLong$IsDMP, ], size = 0.5, color = "grey") +
    geom_vline(data = data.frame(x = regionDMPs$cstart), aes(xintercept = x), color = "lightgrey", linetype = "dashed") +
    geom_rug(data = methylationDataLong[methylationDataLong$IsDMP, ], aes(x = cstart), color = "red", sides = "t") +
    scale_color_manual("", values = c("Control" = "grey", "Case" = "red")) +
    scale_linetype_manual("", values = c("Control" = "solid", "Case" = "solid")) +
    labs(
      title = ICR,
      x = "CpG Coordinates",
      y = "Methylation Value"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  if (interactive) {
    plot <- ggplotly(plot) # Converting to interactive plot
  }
    return(plot)
  
}
