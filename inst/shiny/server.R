# File: inst/shiny/server.R

server <- function(input, output, session) {
  
  # Reactive values to store data
  values <- reactiveValues(
    betaMatrix = NULL,
    sampleInfo = NULL,
    ICRcpg = NULL,
    df.ICR = NULL,
    cpgs_analysis = NULL,
    dmps_results = NULL,
    processed = FALSE
  )
  
  # Check if data is loaded
  output$dataLoaded <- reactive({
    !is.null(values$betaMatrix)
  })
  outputOptions(output, "dataLoaded", suspendWhenHidden = FALSE)
  
  # Check if DMPs analysis is done
  output$dmpsDone <- reactive({
    !is.null(values$dmps_results)
  })
  outputOptions(output, "dmpsDone", suspendWhenHidden = FALSE)
  
  # File upload and preview
  observeEvent(input$betaMatrix, {
    req(input$betaMatrix)
    
    df <- read.csv(input$betaMatrix$datapath, row.names = 1, check.names = FALSE)
    values$betaMatrix <- as.matrix(df)
    
    # Generate sample group UI
    output$sampleGroupUI <- renderUI({
      samples <- colnames(values$betaMatrix)
      
      tagList(
        selectInput("controlSamples", "Select Control Samples:",
                    choices = samples,
                    multiple = TRUE,
                    selected = NULL
        ),
        
        helpText("Samples not selected as controls will be assigned to the Case group.")
      )
    })
    
    # Data summary
    output$dataSummary <- renderPrint({
      cat("Data dimensions:\n")
      cat(sprintf("  CpGs: %d\n", nrow(values$betaMatrix)))
      cat(sprintf("  Samples: %d\n", ncol(values$betaMatrix)))
      cat("\nSample names:\n")
      cat(paste(colnames(values$betaMatrix), collapse = ", "))
    })
    
    # Data preview
    output$dataPreview <- DT::renderDataTable({
      DT::datatable(
        head(values$betaMatrix, 100),
        options = list(
          pageLength = 10,
          scrollX = TRUE
        )
      )
    })
  })
  
  # Process data
  observeEvent(input$processData, {
    req(values$betaMatrix, input$controlSamples)
    
    showModal(modalDialog(
      title = "Processing Data",
      "Extracting CpGs and creating ICR matrices...",
      footer = NULL
    ))
    
    # Create sample info
    samples <- colnames(values$betaMatrix)
    values$sampleInfo <- ifelse(samples %in% input$controlSamples, "Control", "Case")
    
    # Process data
    tryCatch({
      # Extract CpGs
      values$ICRcpg <- make_cpgs(Bmatrix = values$betaMatrix, bedmeth = input$arrayType)
      
      # Create ICR matrix
      values$df.ICR <- make_ICRs(Bmatrix = values$betaMatrix, bedmeth = input$arrayType)
      
      values$processed <- TRUE
      
      removeModal()
      
      showNotification("Data processed successfully!", type = "message")
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error:", e$message), type = "error", duration = 10)
    })
  })
  
  # CpG Coverage Analysis
  observeEvent(input$runCpGAnalysis, {
    req(values$ICRcpg)
    
    showModal(modalDialog(
      title = "Running Analysis",
      "Calculating CpG coverage...",
      footer = NULL
    ))
    
    tryCatch({
      values$cpgs_analysis <- plot_cpgs_coverage(values$ICRcpg, bedmeth = input$arrayType)
      
      removeModal()
      
      # Render plots
      output$cpgCountPlot <- renderPlotly({
        ggplotly(values$cpgs_analysis$plot_counts)
      })
      
      output$cpgPercentagePlot <- renderPlotly({
        ggplotly(values$cpgs_analysis$plot_percentage)
      })
      
      output$cpgCoverageTable <- DT::renderDataTable({
        DT::datatable(
          values$cpgs_analysis$data,
          options = list(pageLength = 25)
        )
      })
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error:", e$message), type = "error", duration = 10)
    })
  })
  
  # ICR Heatmap
  # ICR Heatmap
  observeEvent(input$generateHeatmap, {
    req(values$df.ICR, values$sampleInfo)
    
    # Add validation
    if (is.null(values$df.ICR) || nrow(values$df.ICR) == 0) {
      showNotification("No ICR data available. Please process data first.", type = "error")
      return()
    }
    
    if (is.null(values$sampleInfo) || length(values$sampleInfo) != ncol(values$df.ICR)) {
      showNotification("Sample information doesn't match data dimensions.", type = "error")
      return()
    }
    
    showModal(modalDialog(
      title = "Generating Heatmap",
      "Creating ICR methylation heatmap...",
      footer = NULL
    ))
    
    tryCatch({
      # Generate the heatmap
      p <- ICRs_heatmap(
        df_ICR = values$df.ICR,
        sampleInfo = values$sampleInfo,
        control_label = "Control",
        case_label = "Case",
        bedmeth = input$arrayType,
        order_by = input$orderBy,
        plot_type = input$plotType,
        sd_threshold = ifelse(input$plotType == "defect", input$sdThreshold, 3)
      )
      
      # Render the plot
      output$icrHeatmap <- renderPlot({
        print(p)
      }, height = 800)
      
      removeModal()
      showNotification("Heatmap generated successfully!", type = "message")
      
    }, error = function(e) {
      removeModal()
      
      # More detailed error reporting
      error_msg <- paste("Heatmap generation failed:", e$message)
      showNotification(error_msg, type = "error", duration = 15)
      
      # Print detailed error to console for debugging
      cat("=== HEATMAP ERROR DETAILS ===\n")
      cat("Error message:", e$message, "\n")
      cat("df_ICR dimensions:", dim(values$df.ICR), "\n")
      cat("sampleInfo:", paste(values$sampleInfo, collapse = ", "), "\n")
      cat("arrayType:", input$arrayType, "\n")
      cat("orderBy:", input$orderBy, "\n")
      cat("plotType:", input$plotType, "\n")
      if (input$plotType == "defect") {
        cat("sdThreshold:", input$sdThreshold, "\n")
      }
      cat("============================\n")
      print(traceback())
    })
  })
  
  # Download heatmap - also fix this
  output$downloadHeatmap <- downloadHandler(
    filename = function() {
      paste0("ICR_heatmap_", input$plotType, "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      tryCatch({
        pdf(file, width = 12, height = 10)
        p <- ICRs_heatmap(
          df_ICR = values$df.ICR,
          sampleInfo = values$sampleInfo,
          control_label = "Control",
          case_label = "Case",
          bedmeth = input$arrayType,
          order_by = input$orderBy,
          plot_type = input$plotType,
          sd_threshold = ifelse(input$plotType == "defect", input$sdThreshold, 3)
        )
        print(p)
        dev.off()
      }, error = function(e) {
        showNotification(paste("Download failed:", e$message), type = "error", duration = 10)
      })
    }
  )
  
  # DMP Analysis
  observeEvent(input$runDMPAnalysis, {
    req(values$ICRcpg, values$sampleInfo)
    
    showModal(modalDialog(
      title = "Running DMP Analysis",
      "Identifying differentially methylated positions...",
      footer = NULL
    ))
    
    tryCatch({
      values$dmps_results <- iDMPs(
        data = values$ICRcpg,
        sampleInfo = values$sampleInfo,
        pValueCutoff = input$pValueCutoff
      )
      
      removeModal()
      
      # Update ICR choices for region explorer
      updateSelectInput(session, "selectedICR",
                        choices = unique(values$dmps_results$topDMPs$ICR)
      )
      
      # Render results
      output$dmpsTable <- DT::renderDataTable({
        DT::datatable(
          values$dmps_results$topDMPs,
          options = list(
            pageLength = 25,
            scrollX = TRUE
          )
        ) %>%
          DT::formatRound(columns = c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val"), digits = 4)
      })
      
      # Volcano plot
      output$volcanoPlot <- renderPlotly({
        plot_data <- values$dmps_results$topDMPs
        plot_data$neg_log_p <- -log10(plot_data$P.Value)
        plot_data$significant <- plot_data$P.Value < input$pValueCutoff
        
        p <- ggplot(plot_data, aes(x = logFC, y = neg_log_p, color = significant)) +
          geom_point(alpha = 0.6) +
          scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
          geom_hline(yintercept = -log10(input$pValueCutoff), linetype = "dashed") +
          geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
          labs(x = "Log Fold Change", y = "-log10(p-value)", title = "Volcano Plot") +
          theme_bw()
        
        ggplotly(p)
      })
      
      # Summary
      output$dmpsSummary <- renderPrint({
        cat("DMP Analysis Summary\n")
        cat("==================\n\n")
        cat(sprintf("Total CpGs tested: %d\n", nrow(values$dmps_results$eBayesfit$coefficients)))
        cat(sprintf("Significant DMPs (p < %g): %d\n", input$pValueCutoff, nrow(values$dmps_results$topDMPs)))
        cat(sprintf("\nTop 5 DMPs by p-value:\n"))
        print(head(values$dmps_results$topDMPs[order(values$dmps_results$topDMPs$P.Value), c("ICR", "logFC", "P.Value")], 5))
      })
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error:", e$message), type = "error", duration = 10)
    })
  })
  
  # Download DMP results
  output$downloadDMPs <- downloadHandler(
    filename = function() {
      paste0("DMP_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(values$dmps_results$topDMPs, file, row.names = FALSE)
    }
  )
  
  # Region Explorer
  observeEvent(input$plotRegion, {
    req(input$selectedICR, values$dmps_results, values$ICRcpg)
    
    # Validation checks
    if (is.null(input$selectedICR) || input$selectedICR == "") {
      showNotification("Please select an ICR first.", type = "warning")
      return()
    }
    
    if (nrow(values$dmps_results$topDMPs) == 0) {
      showNotification("No significant DMPs available for plotting.", type = "warning")
      return()
    }
    
    showModal(modalDialog(
      title = "Generating Plot",
      "Creating region plot...",
      footer = NULL
    ))
    
    # Clear previous plots
    output$regionPlot <- renderPlot({ NULL })
    output$regionPlotUI <- renderUI({ NULL })
    
    tryCatch({
      # Debug information
      cat("=== REGION PLOT DEBUG ===\n")
      cat("Selected ICR:", input$selectedICR, "\n")
      cat("Plot interactive:", input$plotInteractive, "\n")
      cat("Available ICRs in DMPs:", paste(unique(values$dmps_results$topDMPs$ICR), collapse = ", "), "\n")
      cat("ICRcpg dimensions:", dim(values$ICRcpg), "\n")
      cat("Sample info length:", length(values$sampleInfo), "\n")
      
      # Check if selected ICR exists in DMPs
      if (!input$selectedICR %in% values$dmps_results$topDMPs$ICR) {
        stop(paste("Selected ICR", input$selectedICR, "not found in significant DMPs"))
      }
      
      # Create the plot
      plot_interactive <- as.logical(input$plotInteractive)
      
      plot <- plot_line_ICR(
        significantDMPs = values$dmps_results$topDMPs,
        ICRcpg = values$ICRcpg,
        ICR = input$selectedICR,
        sampleInfo = values$sampleInfo,
        interactive = plot_interactive
      )
      
      removeModal()
      
      # Render the appropriate UI based on plot type
      if (plot_interactive) {
        output$regionPlotUI <- renderUI({
          tagList(
            h4(paste("Interactive Plot:", input$selectedICR)),
            plotlyOutput("regionPlotInteractive", height = "600px"),
            br(),
            helpText("Use plotly controls to zoom, pan, and hover over data points.")
          )
        })
        
        output$regionPlotInteractive <- renderPlotly({
          plot
        })
      } else {
        output$regionPlotUI <- renderUI({
          tagList(
            h4(paste("Static Plot:", input$selectedICR)),
            plotOutput("regionPlotStatic", height = "600px")
          )
        })
        
        output$regionPlotStatic <- renderPlot({
          print(plot)
        }, height = 600)
      }
      
      showNotification("Region plot generated successfully!", type = "message")
      
    }, error = function(e) {
      removeModal()
      
      error_msg <- paste("Region plot failed:", e$message)
      showNotification(error_msg, type = "error", duration = 15)
      
      # Detailed error logging
      cat("=== REGION PLOT ERROR ===\n")
      cat("Error message:", e$message, "\n")
      cat("Selected ICR:", input$selectedICR, "\n")
      cat("Interactive setting:", input$plotInteractive, "\n")
      cat("DMPs available:", nrow(values$dmps_results$topDMPs), "\n")
      
      if (exists("values$dmps_results$topDMPs")) {
        dmps_for_icr <- values$dmps_results$topDMPs[values$dmps_results$topDMPs$ICR == input$selectedICR, ]
        cat("DMPs for selected ICR:", nrow(dmps_for_icr), "\n")
      }
      
      if (exists("values$ICRcpg")) {
        icr_cpgs <- values$ICRcpg[values$ICRcpg$ICR == input$selectedICR, ]
        cat("CpGs for selected ICR:", nrow(icr_cpgs), "\n")
      }
      
      cat("Stack trace:\n")
      print(traceback())
      cat("========================\n")
      
      # Show fallback message
      output$regionPlotUI <- renderUI({
        div(
          style = "text-align: center; padding: 50px;",
          h4("Plot Generation Failed", style = "color: red;"),
          p("Error:", e$message),
          p("Please try selecting a different ICR or check the console for detailed error information.")
        )
      })
    })
  })
  
  # Add a reset button for clearing plots
  observeEvent(input$selectedICR, {
    # Clear plots when ICR selection changes
    output$regionPlotUI <- renderUI({
      div(
        style = "text-align: center; padding: 50px;",
        h4("Select Plot Type and Click 'Plot Region'"),
        p("Choose an ICR and plot type, then click the 'Plot Region' button to generate the visualization.")
      )
    })
  })
}