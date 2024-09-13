# Function to plot the number of CpGs covered per ICR, comparing df.ICR.cpg with bedmeth
plot_CpG_coverage <- function(df_ICR_cpg, bedmeth_version = "v1") {
  
  # Load the bedmeth data based on the version specified
  if (bedmeth_version == "v1") {
    message("Loading bedEPICv1...")
    data(bedEPICv1)
    bedmeth_data <- bedEPICv1
  } else if (bedmeth_version == "v2") {
    message("Loading bedEPICv2...")
    data(bedEPICv2)
    bedmeth_data <- bedEPICv2
  } else if (bedmeth_version == "450k") {
    message("Loading bed450k...")
    data(bed450k)
    bedmeth_data <- bed450k
  } else {
    stop("Invalid bedmeth version. Choose from 'v1', 'v2', or '450k'.")
  }
  
  # Count the number of CpGs per ICR in df_ICR_cpg (input data)
  df_ICR_cpg_counts <- df_ICR_cpg %>%
    group_by(ICR) %>%
    summarize(n_cpgs = n())
  
   # Load the appropriate ICRs data based on bedmeth input
  if (bedmeth == "v1" || bedmeth == "450k") {
    message("Loading DMRs.hg19...")
    data(DMRs.hg19)
    ICRs <- DMRs.hg19
  } else if (bedmeth == "v2") {
    message("Loading DMRs.hg38...")
    data(DMRs.hg38)
    ICRs <- DMRs.hg38
  }

  bedmeth_counts <- bed_intersect(ICRs, bedmeth_data) %>%
  mutate(chr = gsub("chr", "", chrom)) %>%
  mutate(chr = as.numeric(chr)) %>%
  arrange(chr, start.x) %>%
  dplyr::select(probeID.y, start.y, ICR.x, start.x, end.x) %>%
  as.data.frame() %>%
  dplyr::count(ICR.x, sort = T) %>%
  dplyr::rename("ICR" = "ICR.x",
                "Total_cov" = "n")

  comparison_data <- full_join(df_ICR_cpg_counts, probeICR, by = "ICR")

  comparison_data_long <- comparison_data %>%
  pivot_longer(cols = c(n_cpgs, Total_cov), names_to = "Measure", values_to = "Count") %>%
  replace_na(list(Count = 0))

  ggplot(comparison_data_long, aes(x = reorder(ICR,Count), y = Count, fill = Measure)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5)) +
  labs(title = "Number of CpGs Covered per ICR", x = "ICR", y = "Number of CpGs") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw() + 
  coord_flip()

  return(comparison_data)
