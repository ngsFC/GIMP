#' GIMP Shiny Application
#'
#' @description Launch the GIMP Shiny app for interactive analysis of ICR methylation
#' @export
#' @examples
#' if(interactive()){
#'   GIMP_app()
#' }

GIMP_app <- function() {
  appDir <- system.file("shiny", package = "GIMP")
  if (appDir == "") {
    stop("Could not find Shiny app directory. Try re-installing `GIMP`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}

# File: inst/shiny/app.R
library(shiny)
library(shinydashboard)
library(DT)
library(GIMP)
library(tidyverse)
library(plotly)

# Source UI and server modules
source("ui.R")
source("server.R")

# Run the application
shinyApp(ui = ui, server = server)