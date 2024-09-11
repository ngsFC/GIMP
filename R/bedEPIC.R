# Load required libraries
library(dplyr)
library(tidyr)
library(readr)

# Function to create bedEPIC from EPICv1 or EPICv2
create_bedEPIC <- function(version = "v1") {
  if (version == "v1") {
    # Load EPICv1 annotation
    annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    bedEPIC <- as.data.frame(annEPIC) %>%
      rownames_to_column("probeID") %>%
      dplyr::select(chr, pos, probeID) %>%
      mutate(end = pos) %>%
      relocate(end, .after = pos)

    colnames(bedEPIC) <- c("chrom", "start", "end", "probeID")
    
  } else if (version == "v2") {
    # File URL and destination path
    url <- "https://zwdzwd.github.io/InfiniumAnnotation/EPICv2.hg38.manifest.gencode.v41.tsv.gz"
    destfile <- "data/EPICv2.hg38.manifest.gencode.v41.tsv.gz"
    
    # Check if the file exists locally, if not, download it
    if (!file.exists(destfile)) {
      cat("Downloading EPICv2 manifest file...\n")
      download.file(url, destfile, mode = "wb")
      cat("Download complete.\n")
    }
    
    # Read and unzip the file
    EPICv2_manifest <- read_tsv(destfile) %>%
      as.data.frame()
    
    # Create bedEPIC for EPICv2
    bedEPIC <- as.data.frame(EPICv2_manifest) %>%
      rownames_to_column("probeID") %>%
      dplyr::select(CHR, MAPINFO, probeID) %>%
      mutate(end = MAPINFO) %>%
      dplyr::rename(chrom = CHR, start = MAPINFO) %>%
      relocate(end, .after = start)
  } else {
    stop("Invalid version. Use 'v1' or 'v2'.")
  }
  
  return(bedEPIC)
}

# Example usage
bedEPIC <- create_bedEPIC(version = "v2")
write.csv(bedEPIC, "bedEPIC.csv", row.names = FALSE)

