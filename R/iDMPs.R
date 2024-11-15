iDMPs <- function(data, sampleInfo, pValueCutoff = 0.05) {
  # Load necessary libraries
  library(limma)
  
  # Ensure sampleInfo is a factor
  sampleInfo <- factor(sampleInfo)
  
  # Automatically reorder levels to have "Control" as the reference level
  if (!"Control" %in% levels(sampleInfo)) {
    stop("The 'Control' group must be present in sampleInfo.")
  }
  sampleInfo <- relevel(sampleInfo, ref = "Control")
  
  # Validate and split data into methylation values and additional columns
  if (!is.matrix(data)) {
    stop("The input data must be a matrix or data frame.")
  }
  
  # Extract methylation data (assumed all columns except the last four)
  if (ncol(data) < 5) {
    stop("The data must have at least four additional columns: 'cstart', 'ICR', 'start', 'end'.")
  }
  
  methylationData <- as.matrix(data[, -c((ncol(data) - 3):ncol(data))])
  additionalColumns <- data[, (ncol(data) - 3):ncol(data)]
  
  expectedColumns <- c("cstart", "ICR", "start", "end")
  if (!all(expectedColumns %in% colnames(additionalColumns))) {
    stop(paste("The last four columns of the data must be named:", paste(expectedColumns, collapse = ", ")))
  }
  
  betaToM <- function(beta) {
    if (any(beta <= 0 | beta >= 1)) {
      stop("Beta values must be between 0 and 1 (exclusive).")
    }
    log2(beta / (1 - beta))
  }
  
  mValues <- betaToM(methylationData)
  
  design <- model.matrix(~sampleInfo)
  fit <- lmFit(mValues, design)
  eBayesfit <- eBayes(fit)
  
  groupLabels <- levels(sampleInfo)
  coefName <- paste0("sampleInfo", groupLabels[2])  # Select second group relative to "Control"
  
  if (!coefName %in% colnames(fit$coefficients)) {
    stop(paste("Coefficient", coefName, "not found in the fit model. Check group labels."))
  }
  
  topDMPs <- topTable(eBayesfit, number = Inf, p.value = pValueCutoff, coef = coefName)
  topDMPs <- cbind(topDMPs, additionalColumns[rownames(topDMPs), , drop = FALSE])
  
  return(list(
    fit = fit,
    eBayesfit = eBayesfit,
    topDMPs = topDMPs,
    groupLabels = groupLabels
  ))
}

