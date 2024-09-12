# Load required packages
library(usethis)
library(devtools)
setwd("/media/ciccio/dati/MethArray/GIMP/data-raw")
#usethis::create_project("/media/ciccio/dati/MethArray/ImpDef")

# Create the data directory if it doesn't exist
usethis::use_data_raw()

# Read the data
dmr38 <- read.table("DMRs.hg38.bed", header = TRUE, sep = "\t")
dmr19 <- read.table("DMRs.hg19.bed", header = TRUE, sep = "\t")
bedEPICv1 <- read.table("bedEPICv1.csv", header = TRUE, sep = ",")
bedEPICv2 <- read.table("bedEPICv2.csv", header = TRUE, sep = ",")
bed450k <- read.table("bed450k.csv", header = TRUE, sep = ",")

# Save the data as an .rda file
usethis::use_data(dmr38, overwrite = T)
usethis::use_data(dmr19, overwrite = T)
usethis::use_data(bedEPICv1, overwrite = T)
usethis::use_data(bedEPICv2, overwrite = T)
usethis::use_data(bed450k, overwrite = T)
