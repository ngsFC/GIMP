# GIMP - Genomic Imprinting Methylation Package

[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://cran.r-project.org/)

**GIMP** is an R package designed for analyzing methylation data with a focus on genomic imprinting regions. The package offers tools for handling data from the **EPIC** methylation array (v1 and v2) and provides functions for identifying and analyzing differentially methylated regions (DMRs).

## Features

- Supports both **EPICv1** and **EPICv2** arrays.
- Easy-to-use functions for performing bed file intersections with imprinting control regions (ICRs).
- Automates the process of working with large methylation datasets and helps identify CpGs in specific DMRs.
- Integration with popular bioinformatics tools for methylation analysis.

## Installation

To install the package from GitHub, you will need to use the `devtools` package in R. Follow the steps below to install **GIMP**.

### Step 1: Install the required packages

Make sure you have the necessary R packages installed:

```r
install.packages("devtools")
install.packages("dplyr")
install.packages("tidyr")
install.packages("readr")
```

### Step 2: Install the GIMP package from GitHub

You can install the GIMP package directly from this repository using devtools:

```r
# Load devtools
library(devtools)

# Install the GIMP package from GitHub
install_github("ngsFC/GIMP")

# load the package
library(GIMP)
```

### Usage

1. Loading Imprinting Control Regions (ICRs)

Before analyzing the methylation data, load the DMRs file (ICRs) depending on whether you're working with EPICv1 or EPICv2 data. Here's an example for loading ICRs in your environment:


```r
# For EPICv1
ICRs <- read.table(system.file("data", "DMRs.hg19.bed", package = "GIMP"), header = FALSE)
colnames(ICRs) <- c("chrom", "start", "end", "germ", "ICR")

# For EPICv2
ICRs <- read.table(system.file("data", "DMRs.hg38.bed", package = "GIMP"), header = FALSE)
colnames(ICRs) <- c("chrom", "start", "end", "germ", "ICR")
```

2. Creating the EPIC Bed File

To create the bed file from the EPICv1 or EPICv2 manifest, run the create_bedEPIC function:

```r
# Create bedEPIC for EPICv1
bedEPIC <- create_bedEPIC(version = "v1")

# Create bedEPIC for EPICv2
bedEPIC <- create_bedEPIC(version = "v2")
```

3. Analyzing Methylation Data

After loading the ICRs and creating the bedEPIC, you can use the function create_ICR_matrices() to analyze your methylation data (e.g., myCombat).

```r
# Example of using the function to create ICR matrices
result <- create_ICR_matrices(myCombat = your_combat_data)

# Access the data frames
df.ICR.cpg <- result$df.ICR.cpg
df.ICR <- result$df.ICR
```

### Contributing

We welcome contributions! If you'd like to contribute to this project, please fork the repository and submit a pull request with your changes.

### License

This project is licensed under the MIT License - see the LICENSE file for details.
