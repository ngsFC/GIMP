# Function to load the appropriate DMR file based on EPIC version
load_DMRs <- function(version = "v1") {
  if (version == "v1") {
    # Access DMRs.hg19.bed from the package's inst folder
    file_path <- system.file("data", "DMRs.hg19.bed", package = "GIMP")
    ICRs <- read.table(file_path, header = FALSE)
    colnames(ICRs) <- c("chrom", "start", "end", "germ", "ICR")
  } else if (version == "v2") {
    # Access DMRs.hg38.bed from the package's inst folder
    file_path <- system.file("data", "DMRs.hg38.bed", package = "GIMP")
    ICRs <- read.table(file_path, header = FALSE)
    colnames(ICRs) <- c("chrom", "start", "end", "germ", "ICR")
  } else {
    stop("Invalid version. Use 'v1' or 'v2'.")
  }
  
  # Clean the ICRs dataframe
  ICRs <- ICRs %>%
    mutate(chr = gsub("chr", "", chrom)) %>%
    mutate(chr = as.numeric(chr)) %>%
    arrange(chr, start) %>%
    dplyr::select(-chr)
  
  return(ICRs)
}

