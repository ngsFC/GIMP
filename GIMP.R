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

setwd("/media/ciccio/dati/MethArray/Public_BWS")
df <- readRDS("BWSimp.rds")

ICRcpg <- make_cpgs(Bmatrix = df, bedmeth = "v1")

cpgs_analysis <- plot_CpG_coverage(ICRcpg, bedmeth = "v1")

df.ICR <- make_ICRs(Bmatrix = df, bedmeth = "v1")

icr.heats <- GIMP::iDMR_heatmap(df.ICR, 
                                   group_vector = c(rep("Case", 13), rep("Control", 24)),
                                   control_label = "Control", 
                                   case_label = "Case",
                                   cluster_by = "meth")

