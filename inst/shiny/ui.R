# File: inst/shiny/ui.R

ui <- dashboardPage(
  dashboardHeader(title = "GIMP: Genomic Imprinting Methylation Patterns"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Upload", tabName = "upload", icon = icon("upload")),
      menuItem("CpG Analysis", tabName = "cpg", icon = icon("dna")),
      menuItem("ICR Analysis", tabName = "icr", icon = icon("chart-bar")),
      menuItem("Differential Analysis", tabName = "dmps", icon = icon("chart-line")),
      menuItem("Region Explorer", tabName = "explorer", icon = icon("search")),
      menuItem("Help", tabName = "help", icon = icon("question-circle"))
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper, .right-side {
          background-color: #f4f4f4;
        }
        .box-header {
          background-color: #3c8dbc;
          color: white;
        }
      "))
    ),
    
    tabItems(
      # Data Upload Tab
      tabItem(tabName = "upload",
        fluidRow(
          box(
            title = "Upload Methylation Data",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            
            fluidRow(
              column(6,
                fileInput("betaMatrix", 
                  "Choose Beta Matrix CSV File",
                  accept = c(".csv", ".txt"),
                  placeholder = "Rows: CpGs, Columns: Samples"
                ),
                
                radioButtons("arrayType", "Array Type:",
                  choices = list(
                    "EPIC v1 (hg19)" = "v1",
                    "EPIC v2 (hg38)" = "v2",
                    "450k (hg19)" = "450k"
                  ),
                  selected = "v1"
                ),
                
                helpText("Upload a CSV file where rows are CpG probes and columns are samples. 
                         First column should contain probe IDs.")
              ),
              
              column(6,
                h4("Sample Group Assignment"),
                uiOutput("sampleGroupUI"),
                
                actionButton("processData", "Process Data", 
                  class = "btn-primary",
                  icon = icon("play")
                )
              )
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Data Summary",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            
            verbatimTextOutput("dataSummary"),
            
            DT::dataTableOutput("dataPreview")
          )
        )
      ),
      
      # CpG Analysis Tab
      tabItem(tabName = "cpg",
        fluidRow(
          box(
            title = "CpG Coverage Analysis",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            
            conditionalPanel(
              condition = "output.dataLoaded",
              
              actionButton("runCpGAnalysis", "Generate CpG Analysis", 
                class = "btn-success", 
                icon = icon("calculator")
              ),
              
              br(), br(),
              
              tabsetPanel(
                tabPanel("Coverage Counts",
                  plotlyOutput("cpgCountPlot", height = "600px")
                ),
                tabPanel("Coverage Percentage",
                  plotlyOutput("cpgPercentagePlot", height = "600px")
                ),
                tabPanel("Coverage Data",
                  DT::dataTableOutput("cpgCoverageTable")
                )
              )
            ),
            
            conditionalPanel(
              condition = "!output.dataLoaded",
              h4("Please upload data first in the Data Upload tab.")
            )
          )
        )
      ),
      
      # ICR Analysis Tab
      tabItem(tabName = "icr",
        fluidRow(
          box(
            title = "ICR Methylation Heatmap",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            
            conditionalPanel(
              condition = "output.dataLoaded",
              
              fluidRow(
                column(4,
                  selectInput("plotType", "Plot Type:",
                    choices = list(
                      "Beta Values" = "beta",
                      "Delta Beta (vs Control)" = "delta",
                      "Defect Matrix" = "defect"
                    ),
                    selected = "beta"
                  )
                ),
                
                column(4,
                  selectInput("orderBy", "Order By:",
                    choices = list(
                      "Coordinates" = "cord",
                      "Methylation Values" = "meth"
                    ),
                    selected = "cord"
                  )
                ),
                
                column(4,
                  conditionalPanel(
                    condition = "input.plotType == 'defect'",
                    numericInput("sdThreshold", "SD Threshold:",
                      value = 3, min = 1, max = 5, step = 0.5
                    )
                  )
                )
              ),
              
              actionButton("generateHeatmap", "Generate Heatmap", 
                class = "btn-success", 
                icon = icon("fire")
              ),
              
              br(), br(),
              
              plotOutput("icrHeatmap", height = "800px"),
              
              br(),
              
              downloadButton("downloadHeatmap", "Download Heatmap")
            ),
            
            conditionalPanel(
              condition = "!output.dataLoaded",
              h4("Please upload data first in the Data Upload tab.")
            )
          )
        )
      ),
      
      # Differential Analysis Tab
      tabItem(tabName = "dmps",
        fluidRow(
          box(
            title = "Differential Methylation Analysis",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            
            conditionalPanel(
              condition = "output.dataLoaded",
              
              fluidRow(
                column(6,
                  numericInput("pValueCutoff", "P-value Cutoff:",
                    value = 0.05, min = 0.001, max = 0.1, step = 0.01
                  )
                ),
                
                column(6,
                  actionButton("runDMPAnalysis", "Run DMP Analysis", 
                    class = "btn-success", 
                    icon = icon("calculator")
                  )
                )
              ),
              
              br(),
              
              tabsetPanel(
                tabPanel("Significant DMPs",
                  DT::dataTableOutput("dmpsTable")
                ),
                tabPanel("Volcano Plot",
                  plotlyOutput("volcanoPlot", height = "600px")
                ),
                tabPanel("Summary Statistics",
                  verbatimTextOutput("dmpsSummary")
                )
              ),
              
              br(),
              
              downloadButton("downloadDMPs", "Download DMP Results")
            ),
            
            conditionalPanel(
              condition = "!output.dataLoaded",
              h4("Please upload data first in the Data Upload tab.")
            )
          )
        )
      ),
      
      # Region Explorer Tab
      tabItem(tabName = "explorer",
        fluidRow(
          box(
            title = "ICR Region Explorer",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            
            conditionalPanel(
              condition = "output.dmpsDone",
              
              fluidRow(
                column(6,
                  selectInput("selectedICR", "Select ICR:",
                    choices = NULL
                  )
                ),
                
                column(6,
                  radioButtons("plotInteractive", "Plot Type:",
                    choices = list(
                      "Interactive" = TRUE,
                      "Static" = FALSE
                    ),
                    selected = TRUE
                  )
                )
              ),
              
              actionButton("plotRegion", "Plot Region", 
                class = "btn-success", 
                icon = icon("chart-line")
              ),
              
              br(), br(),
              
              uiOutput("regionPlotUI")
            ),
            
            conditionalPanel(
              condition = "!output.dmpsDone",
              h4("Please run DMP analysis first in the Differential Analysis tab.")
            )
          )
        )
      ),
      
      # Help Tab
      tabItem(tabName = "help",
        fluidRow(
          box(
            title = "GIMP Shiny App Help",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            
            includeMarkdown(system.file("shiny/help.md", package = "GIMP"))
          )
        )
      )
    )
  )
)