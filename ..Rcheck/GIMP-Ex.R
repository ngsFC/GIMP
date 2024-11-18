pkgname <- "GIMP"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('GIMP')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("DMRs.hg19")
### * DMRs.hg19

flush(stderr()); flush(stdout())

### Name: DMRs.hg19
### Title: Imprinted Regions
### Aliases: DMRs.hg19
### Keywords: datasets

### ** Examples

data(DMRs.hg19)
head(DMRs.hg19)



cleanEx()
nameEx("DMRs.hg38")
### * DMRs.hg38

flush(stderr()); flush(stdout())

### Name: DMRs.hg38
### Title: Imprinted Regions
### Aliases: DMRs.hg38
### Keywords: datasets

### ** Examples

data(DMRs.hg38)
head(DMRs.hg38)



cleanEx()
nameEx("bed450k")
### * bed450k

flush(stderr()); flush(stdout())

### Name: bed450k
### Title: BED 450K probes
### Aliases: bed450k
### Keywords: datasets

### ** Examples

data(bed450k)
head(bed450k)



cleanEx()
nameEx("bedEPICv1")
### * bedEPICv1

flush(stderr()); flush(stdout())

### Name: bedEPICv1
### Title: BED EPICv1 probes
### Aliases: bedEPICv1
### Keywords: datasets

### ** Examples

data(bedEPICv1)
head(bedEPICv1)



cleanEx()
nameEx("bedEPICv2")
### * bedEPICv2

flush(stderr()); flush(stdout())

### Name: bedEPICv2
### Title: BED EPICv2 probes
### Aliases: bedEPICv2
### Keywords: datasets

### ** Examples

data(bedEPICv2)
head(bedEPICv2)



cleanEx()
nameEx("create_bedmeth")
### * create_bedmeth

flush(stderr()); flush(stdout())

### Name: create_bedmeth
### Title: Create BED File Data from Methylation Array Annotations
### Aliases: create_bedmeth

### ** Examples

# Create BED-format data with the default version (EPIC v1)
bed_data <- create_bedmeth()
head(bed_data)  # View the first few rows

# Use a different annotation version if available
bed_data_v2 <- create_bedmeth(version = "v2")



cleanEx()
nameEx("iDMR_heatmap")
### * iDMR_heatmap

flush(stderr()); flush(stdout())

### Name: iDMR_heatmap
### Title: Generate Heatmap of Imprinted DMRs Methylation
### Aliases: iDMR_heatmap

### ** Examples

# Example sampleInfo with "Case" and "Control" labels for each sample
sampleInfo <- c(rep("Case", 10), rep("Control", 10))
DMR_heatmap(df_ICR = my_ICR_data, sampleInfo = sampleInfo, annotation_col = list(Sample = c("darkgreen", "darkred")))



cleanEx()
nameEx("make_ICRs")
### * make_ICRs

flush(stderr()); flush(stdout())

### Name: make_ICRs
### Title: Create the ICR Matrix
### Aliases: make_ICRs

### ** Examples

ICRmatrix <- make_ICRs(Bmatrix = df, bedmeth = "v1")



cleanEx()
nameEx("make_cpgs")
### * make_cpgs

flush(stderr()); flush(stdout())

### Name: make_cpgs
### Title: Create ICR CpG Matrix
### Aliases: make_cpgs

### ** Examples

# Generate the ICR CpG matrix with default BED version (EPIC v1)
ICRcpg <- make_cpgs(Bmatrix = df, bedmeth = "v1")

# Use a different BED version, such as EPIC v2
ICRcpg_v2 <- make_cpgs(Bmatrix = df, bedmeth = "v2")



cleanEx()
nameEx("plot_CpG_coverage")
### * plot_CpG_coverage

flush(stderr()); flush(stdout())

### Name: plot_CpG_coverage
### Title: Plot ICR CpG Matrix with Counts and Percentage Coverage
### Aliases: plot_CpG_coverage

### ** Examples

plot_CpG_coverage(df_ICR_cpg_counts, bedmeth = "v1")



cleanEx()
nameEx("plot_line_region")
### * plot_line_region

flush(stderr()); flush(stdout())

### Name: plot_line_region
### Title: Plot Line Plot for Imprinted DMR Methylations
### Aliases: plot_line_region

### ** Examples

# Example data for significantDMPs
plot <- plot_line_region(significantDMPs, ICRcpg, ICR = "KCNQ1OT1:TSS-DMR", sampleInfo = sampleInfo, interactive = T)
print(plot)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
