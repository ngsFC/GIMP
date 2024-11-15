iDMPs <- function(data, sampleInfo, pValueCutoff = 0.001) {
  
  if (!is.factor(sampleInfo)) {
    sampleInfo <- factor(sampleInfo)
  }
  
  sampleInfo <- relevel(sampleInfo, "Control")
  
  betaToM <- function(beta) {
    log2(beta / (1 - beta))
  }
  
  mValues <- betaToM(data)
  design <- model.matrix(~sampleInfo, data = as.data.frame(mValues))
  fit <- lmFit(mValues, design)
  eBayesfit <- eBayes(fit)
  topDMPs <- topTable(eBayesfit, number = Inf, p.value = pValueCutoff, coef = "sampleInfocase")
  
  return(list(
    fit = fit,
    eBayesfit = eBayesfit,
    topDMPs = topDMPs
  ))
}

