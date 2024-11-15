iDMPs <- function(data, sampleInfo, pValueCutoff = 0.05) {
  sampleInfo <- factor(sampleInfo)
  
  if (!"Control" %in% levels(sampleInfo)) {
    stop("The 'Control' group must be present in sampleInfo.")
  }
  sampleInfo <- relevel(sampleInfo, ref = "Control")
  
  if (!is.matrix(data)) {
    stop("The input data must be a matrix.")
  }
  if (!is.numeric(data)) {
    data <- apply(data, 2, as.numeric)
    if (any(is.na(data))) {
      stop("Non-numeric values detected in the data.")
    }
  }
  
  betaToM <- function(beta) {
    if (any(beta <= 0 | beta >= 1)) {
      stop("Beta values must be between 0 and 1 (exclusive).")
    }
    log2(beta / (1 - beta))
  }
  
  mValues <- betaToM(data)
  
  design <- model.matrix(~sampleInfo)
  fit <- lmFit(mValues, design)
  eBayesfit <- eBayes(fit)
  
  groupLabels <- levels(sampleInfo)
  coefName <- paste0("sampleInfo", groupLabels[2])
  
  if (!coefName %in% colnames(fit$coefficients)) {
    stop(paste("Coefficient", coefName, "not found in the fit model. Check group labels."))
  }
  
  topDMPs <- topTable(eBayesfit, number = Inf, p.value = pValueCutoff, coef = coefName)
  
  return(list(
    fit = fit,
    eBayesfit = eBayesfit,
    topDMPs = topDMPs,
    groupLabels = groupLabels
  ))
}

