#' Create the ICR Matrix
#'
#' @examples
#' ICRmatrix <- make_ICRs(Bmatrix = df, bedmeth = "v1")
#' @export

make_ICRs <- function(Bmatrix, bedmeth = "v1") {
  
  # Load the appropriate bedmeth data based on the bedmeth input
  if (bedmeth == "v1") {
    message("Loading bedEPICv1...")
    data(bedEPICv1)
    bedmeth_data <- bedEPICv1
  } else if (bedmeth == "v2") {
    message("Loading bedEPICv2...")
    data(bedEPICv2)
    bedmeth_data <- bedEPICv2
  } else if (bedmeth == "450k") {
    message("Loading bed450k...")
    data(bed450k)
    bedmeth_data <- bed450k
  } else {
    stop("Invalid bedmeth input. Choose from 'v1', 'v2', or '450k'.")
  }
  
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

  # Perform bed intersection between ICRs and bedmeth_data
  probeICR <- bed_intersect(ICRs, bedmeth_data) %>%
    mutate(chr = gsub("chr", "", chrom)) %>%
    mutate(chr = as.numeric(chr)) %>%
    arrange(chr, start.x) %>%
    dplyr::select(probeID.y, start.y, ICR.x, start.x, end.x) %>%
    as.data.frame()

  colnames(probeICR) <- c("probeID", "cstart", "ICR", "start", "end")

  # Create df.ICR matrix
  df.ICR <- as.data.frame(Bmatrix) %>%
    rownames_to_column("probeID") %>%
    full_join(probeICR, by = "probeID") %>%
    filter(!is.na(ICR)) %>%
    filter(!is.na(.[[2]])) %>%  # Dynamically filter based on the first data column
    dplyr::select(-probeID) %>%
    group_by(ICR) %>%
    summarise_all(mean) %>%
    ungroup() %>%
    full_join(ICRs, by = "ICR") %>%
    mutate(chr = gsub("chr", "", chrom)) %>%
    mutate(chr = as.numeric(chr)) %>%
    arrange(chr, start) %>%
    column_to_rownames("ICR") %>%
    filter(!is.na(.[[1]])) %>%  # Dynamically filter based on the first column of the grouped data
    dplyr::select(-c("chrom", "start", "end", "germ", "chr")) %>%
    as.data.frame()

  return(df.ICR)
}
