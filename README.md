# GIMP: Genomic Imprinting Methylation Patterns

[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://cran.r-project.org/)

<div align="center">
  <img src="GIMP.log.png" alt="Logo" width="200"/>
</div>


**GIMP** (Genomic Imprinting Methylation Patterns) is an R package designed for the analysis of ICRs (Imprinting Control Regions) from methylation array data. It provides a pipeline for extract imprinted CpGs (iCpGs), computing coverage, and analyzing ICRs in probe and sample specific manner. The package supports multiple platforms, including Illumina's 450k, EPIC v1, and EPIC v2 arrays.

### Dependencies

The GIMP package depends on the following R packages:

    tidyverse
    valr
    reshape2
    ggplot2
    pheatmap
    viridisLite
    ggplotify
    readr
    shiny
    plotly
    BiocManager
    IlluminaHumanMethylation450kanno.ilmn12.hg19
    IlluminaHumanMethylationEPICanno.ilm10b4.hg19
    IlluminaHumanMethylationEPICv2anno.20a1.hg38
    limma

## Usage

Below is a sample workflow using GIMP to analyze a Betavalue matrix. This example demonstrates how to generate CpG sites, visualize CpG coverage, create ICRs, and analyze ICR regions.

### Load the Betavalue Matrix

Load the Betavalue matrix. The Betavalue matrix should contain methylation data where rows represent CpG sites and columns represent samples.

```r
# Load the Betavalue matrix (replace with your file)
df <- readRDS("example.rds")
```

### Generate iCpG Sites (make_cpgs())

The make_cpgs() function processes the Betavalue matrix to extract the iCpG for further analysis. You need to specify the array version you are using (e.g., "v1" for EPIC v1).

```r
# Generate CpG sites using the specified BedMeth version (e.g., "v1" for EPICv1)
ICRcpg <- make_cpgs(Bmatrix = df, bedmeth = "v1")
```

### Visualize CpG Coverage (plot_cpg_coverage())

Once the iCpG are extracted, use the plot_cpg_coverage() function to visualize the iCpG coverage. This function provides both the coverage data and a plot for interpretation.

```r
# Visualize CpG coverage
cpgs_analysis <- plot_cpgs_coverage(ICRcpg, bedmeth = "v1")

# The result includes both a plot and data, which can be accessed as:
# cpgs_analysis$plot_counts  # to view the probe coverage for each ICRs
# cpgs_analysis$plot_counts  # to view the percentage of probes covered at ICRs
# cpgs_analysis$data  # to view the coverage data
```

### Look at ICRs (make_ICRs())

The make_ICRs() function creates ICR regions from the Betavalue matrix.

```r
# Generate ICRs
df.ICR <- make_ICRs(Bmatrix = df, bedmeth = "v1")
```

### Analyze ICRs

## Heatmaps

The ICRs_heatmap() function generates heatmapa to visualize the methylation patterns of ICRs using beta values, delta-beta normalization, or a defect matrix based on standard deviations. It is designed to handle genomic methylation data with options for customization in row ordering, annotation, and visualization styles. This function is highly versatile, catering to exploratory analysis of epigenomic data across control and case groups.

```r
sampleInfo <- c(rep("Control", 5), rep("Case", 5))

ICRs_heatmap(
  df_ICR = df_ICR,
  sampleInfo = sampleInfo,
  bedmeth = "v1",
  order_by = "meth",
  plot_type = "beta",
  annotation_col = list(Sample = c(Control = "blue", Case = "red"))
)
```

## iDMPs

The iDMPs() function identifies imprinted differentially methylated probes (iDMPs) by leveraging linear modeling and empirical Bayes methods. This function is designed for comparative analysis of methylation levels between a Control group and another experimental group. It outputs the top-ranked DMPs based on a specified p-value cutoff and includes additional genomic metadata for each position.

```r
result <- iDMPs(data = example_data, sampleInfo = sampleInfo, pValueCutoff = 0.05)

# View top DMPs
head(result$topDMPs)
```

## Plot specific ICR

The plot_line_ICR() function creates a line plot to visualize methylation values across a specific ICR. The plot highlights iDMPs and supports both static and interactive visualizations. This tool is ideal for examining methylation patterns across iCpG within a specific genomic region, with group-based color differentiation.

```r
plot <- plot_line_ICR(
  significantDMPs = significantDMPs,
  ICRcpg = ICRcpg,
  ICR = "ICR1", # use an ICRs name
  sampleInfo = sampleInfo,
  interactive = TRUE
)

# Display the plot
plot
```

## Acknowledgments

The GIMP package was developed by Francesco Cecere [ngsFC]. We gratefully acknowledge the "National Centre for HPC, Big Data and Quantum Computing".

## Contact

For questions, bug reports, or feature requests, please open an issue on GitHub or contact [francesco.cecerengs@gmail.com].
