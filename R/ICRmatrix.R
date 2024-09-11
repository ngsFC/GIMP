  #' Create the Imprinted Matrix
  #'
  #' @examples
  #' plot_heat(score_data, "plots/heat_percentiles.png")
  #' @import tidyverse
  #' @import valr
  #' @export

# Define the function to create df.ICR.cpg and df.ICR
create_ICR_matrices <- function(myCombat, bedEPIC_file = "bedEPIC.csv") {
  
  # Load bedEPIC.csv (this should be the file created beforehand)
  bedEPIC <- read.csv(bedEPIC_file)  # Load the previously created bedEPIC file
  
  # Check if ICRs is loaded in the environment
  if (!exists("ICRs")) {
    stop("ICRs data not found in the environment. Please load the appropriate ICRs data before running this function.")
  }

  # Perform bed intersection between ICRs and bedEPIC
  probeICR <- bed_intersect(ICRs, bedEPIC) %>%
    mutate(chr = gsub("chr", "", chrom)) %>%
    mutate(chr = as.numeric(chr)) %>%
    arrange(chr, start.x) %>%
    dplyr::select(probeID.y, start.y, ICR.x, start.x, end.x) %>%
    as.data.frame()

  colnames(probeICR) <- c("probeID", "cstart", "ICR", "start", "end")

  # Create df.ICR.cpg matrix
  df.ICR.cpg <- as.data.frame(myCombat) %>%
    rownames_to_column("probeID") %>%
    full_join(probeICR, by = "probeID") %>%
    na.omit() %>%
    group_by(ICR) %>%
    column_to_rownames("probeID") %>%
    na.omit() %>%
    as.data.frame()

  # Save df.ICR.cpg
  write.csv(df.ICR.cpg, "df_ICR_cpg.csv", row.names = TRUE)
  message("df_ICR_cpg.csv has been saved.")

  # Create df.ICR matrix
  df.ICR <- as.data.frame(myCombat) %>%
    rownames_to_column("probeID") %>%
    full_join(probeICR, by = "probeID") %>%
    filter(!is.na(ICR), !is.na(P1)) %>%
    dplyr::select(-probeID) %>%
    group_by(ICR) %>%
    summarise_all(mean) %>%
    ungroup() %>%
    full_join(ICRs, by = "ICR") %>%
    mutate(chr = gsub("chr", "", chrom)) %>%
    mutate(chr = as.numeric(chr)) %>%
    arrange(chr, start) %>%
    column_to_rownames("ICR") %>%
    filter(!is.na(P1)) %>%
    dplyr::select(-c("chrom", "start", "end", "germ", "chr")) %>%
    as.data.frame()

  # Save df.ICR
  write.csv(df.ICR, "df_ICR.csv", row.names = TRUE)
  message("df_ICR.csv has been saved.")
  
  return(list(df.ICR.cpg = df.ICR.cpg, df.ICR = df.ICR))
}

# Example usage:
# Assuming ICRs and bedEPIC have been loaded and myCombat is provided
# result <- create_ICR_matrices(myCombat = your_combat_data)

