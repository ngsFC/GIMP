# GIMP: Genomic Imprinting Methylation Patterns

[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://cran.r-project.org/)

<div align="center">
  <img src="GIMP.logo.png" alt="Logo" width="200"/>
</div>


**GIMP** (Genomic Imprinting Methylation Patterns) is an R package designed for the analysis of methylation arrays. It provides a pipeline for generating CpG sites, computing CpG coverage, and analyzing ICRs (Imprinting Control Regions) from methylation data. The package supports multiple platforms, including Illumina's 450k, EPIC v1, and EPIC v2 arrays.

## Installation

You can install the latest version of GIMP directly from GitHub using the devtools package:

```r
# Install devtools if you don't have it
install.packages("devtools")

# Install GIMP from GitHub
devtools::install_github("ngsFC/GIMP", force = TRUE)

# Load the GIMP package
library(GIMP)
```

### Dependencies

The GIMP package depends on the following R packages:

    tidyverse: A collection of R packages for data manipulation.
    valr: Tools to work with genome intervals.
    reshape2: A package for reshaping data between wide and long formats.

Ensure that these dependencies are installed before running GIMP.

```r
# Install required packages
install.packages(c("tidyverse", "valr", "reshape2"))
```

## Usage

Below is a sample workflow using GIMP to analyze a Betavalue matrix. This example demonstrates how to generate CpG sites, visualize CpG coverage, create ICRs, and analyze ICR regions.

### Step 1: Load the Betavalue Matrix

Load the Betavalue matrix. The Betavalue matrix should contain methylation data where rows represent CpG sites and columns represent samples.

```r
# Load the Betavalue matrix (replace with your file)
df <- readRDS("example.rds")
```

### Step 2: Generate CpG Sites (make_cpgs())

The make_cpgs() function processes the Betavalue matrix to generate CpG sites for further analysis. You need to specify the array version you are using (e.g., "v1" for EPIC v1).

```r
# Generate CpG sites using the specified BedMeth version (e.g., "v1" for EPICv1)
ICRcpg <- make_cpgs(Bmatrix = df, bedmeth = "v1")
```

### Step 3: Visualize CpG Coverage (plot_CpG_coverage())

Once the CpG sites are generated, use the plot_CpG_coverage() function to visualize the CpG coverage. This function provides both the coverage data and a plot for interpretation.

```r
# Visualize CpG coverage
cpgs_analysis <- plot_CpG_coverage(ICRcpg, bedmeth = "v1")

# The result includes both a plot and data, which can be accessed as:
# cpgs_analysis$plot  # to view the plot
# cpgs_analysis$data  # to view the coverage data
```

### Step 4: Create Imprinting Control Regions (ICRs) (make_ICRs())

The make_ICRs() function creates ICR regions from the Betavalue matrix.

```r
# Generate ICRs
df.ICR <- make_ICRs(Bmatrix = df, bedmeth = "v1")
```

### Step 5: Analyze ICRs (analyze_ICR())

The analyze_ICR() function performs a differential analysis on the ICRs. You will need to provide a vector that indicates the group assignment (e.g., "Case" or "Control") for each sample.

```r
# Perform ICR analysis
group_vector <- c(rep("Case", 13), rep("Control", 24))  # Define cases and controls
icr_analysis <- GIMP::analyze_ICR(df.ICR, 
                                  group_vector = group_vector, 
                                  control_label = "Control", 
                                  case_label = "Case")
```

### Step 6: Output and Interpretation

After running the ICR analysis, analyze_ICR() returns three heatmaps to help visualize the differences between cases and controls in terms of ICR methylation patterns.

```r
# Access the results (heatmaps, data summaries, etc.)
icr_analysis$heatmap1  # First heatmap
icr_analysis$heatmap2  # Second heatmap
icr_analysis$heatmap3  # Third heatmap
```

## Functions Overview

1. make_cpgs(Bmatrix, bedmeth)

    Description: Generates CpG sites from the Betavalue matrix.
    Arguments:
        Bmatrix: Betavalue matrix (CpG sites by samples).
        bedmeth: Version of the array platform ("v1" for EPICv1, "v2" for EPICv2, etc.).

2. plot_CpG_coverage(cpgs_data, bedmeth)

    Description: Visualizes CpG coverage and returns a plot and data.
    Arguments:
        cpgs_data: Data generated from make_cpgs().
        bedmeth: Version of the array platform.

3. make_ICRs(Bmatrix, bedmeth)

    Description: Creates ICRs from the Betavalue matrix.
    Arguments:
        Bmatrix: Betavalue matrix (CpG sites by samples).
        bedmeth: Version of the array platform.

4. analyze_ICR(ICRmatrix, group_vector, control_label, case_label)

    Description: Performs differential methylation analysis on ICRs and returns heatmaps.
    Arguments:
        ICRmatrix: ICR matrix.
        group_vector: Vector indicating the group assignment for each sample (e.g., "Case" or "Control").
        control_label: Label for the control group.
        case_label: Label for the case group.

## License

GIMP is released under the MIT License.

## Acknowledgments

The GIMP package was developed by Francesco Cecere [ngsFC] as part of a project on genomic imprinting and methylation analysis. We gratefully acknowledge the contributions of the R and bioinformatics communities for the packages used in this project.

## Contact

For questions, bug reports, or feature requests, please open an issue on GitHub or contact [francesco.cecerengs@gmail.com].
