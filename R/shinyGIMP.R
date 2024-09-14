library(shiny)
library(GIMP)

# Define UI for application
ui <- fluidPage(
    # Application title
    titlePanel("Genomic Data Analysis Shiny App"),
    
    # Sidebar layout for user inputs and action buttons
    sidebarLayout(
        sidebarPanel(
            # Input: File upload for Betavalue matrix
            fileInput("file1", "Upload Betavalue Matrix (.csv or .txt)", 
                      accept = c(".csv", ".txt")),
            
            # Button to run make_cpgs() and make_ICRs()
            actionButton("runICR", "Generate ICR Matrix"),
            
            # Button to run plot_CpG_coverage()
            actionButton("runCpGCoverage", "Plot CpG Coverage"),
            
            # Button to run analyze_ICR()
            actionButton("runAnalyzeICR", "Analyze ICR and Generate Heatmaps"),
            
            # Download button for ICR matrix
            downloadButton("downloadICR", "Download ICR Matrix")
        ),
        
        # Main panel for displaying outputs
        mainPanel(
            h3("Generated ICR Matrix"),
            tableOutput("icrTable"),    # Display the ICR matrix
            plotOutput("cpgCoveragePlot"),  # Display CpG coverage plot
            tableOutput("cpgCoverageDF"),   # Display CpG coverage data frame
            plotOutput("heatmap1"),     # Heatmap 1 from analyze_ICR()
            plotOutput("heatmap2"),     # Heatmap 2 from analyze_ICR()
            plotOutput("heatmap3")      # Heatmap 3 from analyze_ICR()
        )
    )
)

# Define server logic
server <- function(input, output, session) {
    
    # Reactive value to store the uploaded data and ICR matrix
    betavalueData <- reactiveVal()
    icrMatrix <- reactiveVal()

    # Upload the Betavalue matrix
    observeEvent(input$file1, {
        req(input$file1)
        betavalueData(read.csv(input$file1$datapath))
    })
    
    # Generate ICR Matrix using make_cpgs() and make_ICRs()
    observeEvent(input$runICR, {
        req(betavalueData())
        cpgs <- make_cpgs(betavalueData())
        icr_result <- make_ICRs(cpgs)
        icrMatrix(icr_result)
        
        # Show the ICR matrix in the UI
        output$icrTable <- renderTable({
            icrMatrix()
        })
    })
    
    # Download ICR Matrix
    output$downloadICR <- downloadHandler(
        filename = function() {
            paste("ICR_Matrix.csv")
        },
        content = function(file) {
            write.csv(icrMatrix(), file)
        }
    )
    
    # Plot CpG coverage
    observeEvent(input$runCpGCoverage, {
        req(betavalueData())
        cpgCoverage <- plot_CpG_coverage(betavalueData())
        
        # Show the coverage plot and data frame
        output$cpgCoveragePlot <- renderPlot({
            cpgCoverage$plot  # Assuming the function returns a list with plot
        })
        
        output$cpgCoverageDF <- renderTable({
            cpgCoverage$data  # Assuming the function returns a list with data
        })
    })
    
    # Run analyze_ICR and generate heatmaps
    observeEvent(input$runAnalyzeICR, {
        req(icrMatrix())
        analysis_result <- analyze_ICR(icrMatrix())
        
        # Assuming the function returns a list of heatmaps
        output$heatmap1 <- renderPlot({
            analysis_result$heatmap1  # Heatmap 1
        })
        
        output$heatmap2 <- renderPlot({
            analysis_result$heatmap2  # Heatmap 2
        })
        
        output$heatmap3 <- renderPlot({
            analysis_result$heatmap3  # Heatmap 3
        })
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

