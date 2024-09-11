.onLoad <- function(libname, pkgname) {
  # Load necessary libraries
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(tidyr))
  suppressPackageStartupMessages(library(tibble))
  suppressPackageStartupMessages(library(readr))
  suppressPackageStartupMessages(library(valr))
  suppressPackageStartupMessages(library(tidyverse))

  message("GIMP package loaded with necessary libraries.")
}
