#install.packages("/media/ciccio/dati/MethArray/GIMP.tar", repos = NULL, type = "source")
#install.packages("remotes")
remove.packages("GIMP")
rm(list = ls())
remotes::install_github("ngsFC/GIMP", auth_token = "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")

library(GIMP)
library(tidyverse)
library(valr)
library(reshape2)
library(grid)
library(pheatmap)
library(viridisLite)
library(ggplotify)
library(limma)
library(plotly)

setwd("/home/ciccio/Desktop/project/data")
df <- readRDS("../BWSimp.rds")

ICRcpg <- make_cpgs(Bmatrix = df, bedmeth = "v1")

cpgs_analysis <- plot_CpG_coverage(ICRcpg, bedmeth = "v1")
cpgs_analysis$plot_counts
cpgs_analysis$plot_percentage

df.ICR <- make_ICRs(Bmatrix = df, bedmeth = "v1")

iDMR_heatmap(df.ICR, 
             sampleInfo = c(rep("Case", 13), rep("Control", 24)),
             control_label = "Control", 
             case_label = "Case",
             order_by = "meth",
             plot_type = "defect",
             sd_threshold = 3)

dmps <- iDMPs(
  data = ICRcpg,
  sampleInfo = c(rep("Case", 13), rep("Control", 24)), 
  pValueCutoff = 0.05
)

significantDMPs <- dmps$topDMPs
print(significantDMPs)

sampleInfo <- c(rep("Case", 13), rep("Control", 24))

plot <- plot_line_region(significantDMPs, ICRcpg, ICR = "KCNQ1OT1:TSS-DMR", sampleInfo = sampleInfo, interactive = T)
print(plot)
