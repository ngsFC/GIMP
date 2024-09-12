library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(tidyverse)

create_bedmeth <- function(version = "v1") {
  if (version == "v1") {
    anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    bedmeth <- as.data.frame(anno) %>%
      rownames_to_column("probeID") %>%
      dplyr::select(chr, pos, probeID) %>%
      mutate(end = pos) %>%
      relocate(end, .after = pos)
    
    colnames(bedmeth) <- c("chrom", "start", "end", "probeID")
    
  } else if (version == "v2") {
    anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    bedmeth <- as.data.frame(anno) %>%
      rownames_to_column("probeID") %>%
      dplyr::select(chr, pos, probeID) %>%
      mutate(end = pos) %>%
      relocate(end, .after = pos)
    
    colnames(bedmeth) <- c("chrom", "start", "end", "probeID")
  } else if (version == "450k") {
    anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    bedmeth <- as.data.frame(anno) %>%
      rownames_to_column("probeID") %>%
      dplyr::select(chr, pos, probeID) %>%
      mutate(end = pos) %>%
      relocate(end, .after = pos)
    
    colnames(bedmeth) <- c("chrom", "start", "end", "probeID")
    } else {
    stop("Invalid version. Use '450k','v1' or 'v2'.")
  }
  
  return(bedmeth)
  }

bed450k <- create_bedmeth(version = "450k")
write.csv(bed450k, "bed450k.csv", row.names = FALSE)
  
bedEPICv1 <- create_bedmeth(version = "v1")
write.csv(bedEPICv1, "bedEPICv1.csv", row.names = FALSE)

bedEPICv2 <- create_bedmeth(version = "v2")
write.csv(bedEPICv2, "bedEPICv2.csv", row.names = FALSE)

