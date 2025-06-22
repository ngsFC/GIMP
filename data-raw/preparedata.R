# Load required packages
library(usethis)
library(devtools)
library(tidyverse)
setwd("/media/ciccio/dati/MethArray/GIMP/data-raw")

# Create the data directory if it doesn't exist
usethis::use_data_raw()

# Read the data
DMRs.hg38 <- read.table("DMRs.hg38.bed", header = F, sep = "\t")
colnames(DMRs.hg38) <- c("chrom", "start", "end", "ICR")
DMRs.hg19 <- read.table("DMRs.hg19.bed", header = F, sep = "\t") %>% dplyr::select(-V4)
colnames(DMRs.hg19) <- c("chrom", "start", "end", "ICR")
bedEPICv1 <- read.table("bedEPICv1.csv", header = TRUE, sep = ",")
bedEPICv2 <- read.table("bedEPICv2.csv", header = TRUE, sep = ",")
bed450k <- read.table("bed450k.csv", header = TRUE, sep = ",")

# Save the data as an .rda file
usethis::use_data(DMRs.hg38, overwrite = T)
usethis::use_data(DMRs.hg19, overwrite = T)
usethis::use_data(bedEPICv1, overwrite = T)
usethis::use_data(bedEPICv2, overwrite = T)
usethis::use_data(bed450k, overwrite = T)
