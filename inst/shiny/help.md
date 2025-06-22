# GIMP Shiny App User Guide

## Overview
The GIMP (Genomic Imprinting Methylation Patterns) Shiny app provides an interactive interface for analyzing methylation patterns at Imprinted Control Regions (ICRs) from methylation array data.

## Workflow

### 1. Data Upload
- **File Format**: Upload a CSV file containing methylation beta values
  - Rows: CpG probes (with probe IDs as row names)
  - Columns: Samples
  - Values: Beta values (0-1)
  
- **Array Type**: Select your methylation array platform:
  - EPIC v1 (hg19)
  - EPIC v2 (hg38)
  - 450k (hg19)

- **Sample Groups**: Assign samples to Control or Case groups
  - Select control samples from the dropdown
  - Unselected samples will be assigned to the Case group

### 2. CpG Analysis
Analyzes the coverage of CpG probes within ICRs:
- **Coverage Counts**: Bar plot showing probe counts per ICR
- **Coverage Percentage**: Percentage of possible probes covered
- **Coverage Data**: Detailed table with coverage statistics

### 3. ICR Analysis
Generates heatmaps of methylation patterns across all ICRs:
- **Plot Types**:
  - *Beta Values*: Raw methylation values (0-1)
  - *Delta Beta*: Difference from control mean
  - *Defect Matrix*: Binary visualization of abnormal methylation
  
- **Ordering Options**:
  - *Coordinates*: ICRs ordered by genomic position
  - *Methylation*: ICRs clustered by methylation patterns

### 4. Differential Analysis
Identifies differentially methylated positions (DMPs) between groups:
- **P-value Cutoff**: Adjust significance threshold
- **Outputs**:
  - Table of significant DMPs
  - Volcano plot visualization
  - Summary statistics

### 5. Region Explorer
Detailed visualization of specific ICR regions:
- Select an ICR from significant DMPs
- View methylation profiles across samples
- Highlights significant DMPs within the region
- Choose between interactive (plotly) or static plots

## Input File Requirements

### Beta Matrix CSV Format
```
,Sample1,Sample2,Sample3,Sample4
cg00000029,0.461,0.472,0.891,0.825
cg00000165,0.123,0.145,0.178,0.189
cg00000236,0.891,0.872,0.461,0.472
```

### Important Notes
- First column must contain probe IDs
- Values must be between 0 and 1
- No missing values in the beta matrix
- Sample names should not contain special characters

## Interpretation Guide

### Heatmap Patterns
- **Normal Imprinting**: Consistent methylation levels within each group
- **Loss of Imprinting**: Abnormal methylation patterns in cases
- **Mosaic Patterns**: Mixed methylation indicating cellular heterogeneity

### DMP Analysis
- **logFC**: Log fold change in methylation (M-values)
- **P.Value**: Statistical significance of difference
- **adj.P.Val**: Multiple testing corrected p-value

### Common Issues
1. **No ICRs Found**: Ensure probe IDs match the selected array type
2. **Analysis Errors**: Check for missing values or non-numeric data
3. **Memory Issues**: Large datasets may require increased memory allocation

## Export Options
- Heatmaps: PDF format
- DMP Results: CSV format
- Plots: Interactive plots can be saved using plotly tools

## Contact
For bugs or feature requests, please visit the [GIMP GitHub repository](https://github.com/ngsFC/GIMP).