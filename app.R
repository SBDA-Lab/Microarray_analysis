# Microarray Analysis Pipeline - Shiny Interface

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install required Bioconductor packages
required_packages <- c("shinydashboard","DT","plotly")

for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE))
    BiocManager::install(package)
}
library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
source("microarray_pipeline.R")


options(shiny.maxRequestSize = 30*1024^2)
# UI Definition
ui <- dashboardPage(
  skin = "blue",
  
  # Header
  dashboardHeader(title = "Microarray Analysis Pipeline"),
  
  # Sidebar
  dashboardSidebar(
    sidebarMenu(
      menuItem("Upload Data", tabName = "upload", icon = icon("upload")),
      menuItem("Quality Control", tabName = "qc", icon = icon("chart-bar")),
      menuItem("Analysis Results", tabName = "results", icon = icon("table")),
      menuItem("Help", tabName = "help", icon = icon("question-circle"))
    )
  ),
  
  # Main Panel
  dashboardBody(
    tabItems(
      # Upload Data Tab
      tabItem(
        tabName = "upload",
        fluidRow(
          box(
            title = "Data Upload",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            
            # File upload
            fileInput("cel_files", "Upload CEL Files",
                     multiple = TRUE,
                     accept = c(".cel", ".CEL")),
            
            # Group assignment
            uiOutput("group_assignments"),
            
            # Output directory
            textInput("output_dir", "Output Directory", 
                     value = "microarray_results"),
            
            # Run analysis button
            actionButton("run_analysis", "Run Analysis",
                        class = "btn-success")
          )
        ),
        
        # Progress and status messages
        fluidRow(
          box(
            width = 12,
            status = "info",
            uiOutput("status_messages")
          )
        )
      ),
      
      # Quality Control Tab
      tabItem(
        tabName = "qc",
        fluidRow(
          box(
            title = "Quality Control Plots",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            
            tabsetPanel(
              tabPanel("RNA Degradation",
                      plotOutput("rna_deg_plot")),
              tabPanel("Box Plots",
                      plotOutput("boxplot")),
              tabPanel("MA Plots",
                      plotOutput("ma_plot")),
              tabPanel("Density Plots",
                      plotOutput("density_plot"))
            )
          )
        )
      ),
      
      # Results Tab
      tabItem(
        tabName = "results",
        fluidRow(
          box(
            title = "Differential Expression Results",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            
            # Comparison selector
            selectInput("comparison_select", "Select Comparison",
                       choices = NULL),
            
            # Results table
            DTOutput("de_table"),
            
            # Download button
            downloadButton("download_results", "Download Results")
          )
        ),
        fluidRow(
          box(
            title = "Volcano Plot",
            width = 6,
            status = "primary",
            solidHeader = TRUE,
            plotlyOutput("volcano_plot")
          ),
          box(
            title = "Top 25 DE Genes Heatmap",
            width = 6,
            status = "primary",
            solidHeader = TRUE,
            plotOutput("heatmap_plot", height = "600px")
          )
        ),
        fluidRow(
          box(
            title = "Normalized Expression Boxplot",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            plotOutput("norm_boxplot")
          )
        )
      ),
      
      # Help Tab
      tabItem(
        tabName = "help",
        fluidRow(
          box(
            title = "User Guide",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            
            HTML("
              <h4>Getting Started</h4>
              <ol>
                <li>Upload your CEL files using the 'Upload Data' tab</li>
                <li>Assign groups to your samples</li>
                <li>Specify an output directory</li>
                <li>Click 'Run Analysis' to start the pipeline</li>
              </ol>
              
              <h4>Quality Control</h4>
              <p>The QC tab shows various diagnostic plots to assess data quality:</p>
              <ul>
                <li>RNA Degradation Plot: Assess RNA quality</li>
                <li>Box Plots: Compare distribution of expression values</li>
                <li>MA Plots: Assess intensity-dependent bias</li>
                <li>Density Plots: View expression value distributions</li>
              </ul>
              
              <h4>Results</h4>
              <p>The Results tab shows:</p>
              <ul>
                <li>Differential expression analysis results</li>
                <li>Interactive volcano plots</li>
                <li>Options to download results</li>
              </ul>
            ")
          )
        )
      )
    )
  )
)

# Server Definition
server <- function(input, output, session) {
  # Reactive values
  values <- reactiveValues(
    cel_files = NULL,
    results_ready = FALSE,
    current_results = NULL,
    comparisons = NULL
  )
  
  # Dynamic group assignment UI
  output$group_assignments <- renderUI({
    req(input$cel_files)
    file_names <- input$cel_files$name
    
    tagList(
      h4("Assign Groups"),
      lapply(seq_along(file_names), function(i) {
        fluidRow(
          column(6, strong(file_names[i])),
          column(6, 
                 selectInput(paste0("group_", i),
                            label = NULL,
                            choices = c("control", "treatment")))
        )
      })
    )
  })
  
  # Run analysis observer
  observeEvent(input$run_analysis, {
    req(input$cel_files)
    
    # Create progress bar
    progress <- shiny::Progress$new()
    progress$set(message = "Running analysis...", value = 0)
    
    # Get group assignments
    groups <- sapply(seq_along(input$cel_files$name), function(i) {
      input[[paste0("group_", i)]]
    })
    
    # Create temporary directory for CEL files
    temp_dir <- tempdir()
    sapply(seq_along(input$cel_files$datapath), function(i) {
      file.copy(input$cel_files$datapath[i],
                file.path(temp_dir, input$cel_files$name[i]))
    })
    
    # Run pipeline
    tryCatch({
      run_microarray_pipeline(temp_dir, input$output_dir, groups)
      values$results_ready <- TRUE
      
      # Update comparison choices
      values$comparisons <- list.files(file.path(input$output_dir, "DE_Analysis"),
                                     pattern = "DE_results_.*\\.csv",
                                     full.names = TRUE)
      updateSelectInput(session, "comparison_select",
                       choices = basename(values$comparisons))
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message),
                      type = "error")
    })
    
    progress$close()
  })
  
  # Load and display results
  observe({
    req(values$results_ready, input$comparison_select)
    
    results_file <- file.path(input$output_dir, "DE_Analysis",
                             input$comparison_select)
    values$current_results <- read.csv(results_file)
    
    # Update results table
    output$de_table <- renderDT({
      datatable(values$current_results,
                options = list(pageLength = 10,
                             scrollX = TRUE))
    })
    
    # Update volcano plot
    output$volcano_plot <- renderPlotly({
      plot_ly(data = values$current_results,
              x = ~logFC,
              y = ~-log10(adj.P.Val),
              text = ~Symbol,
              mode = "markers",
              marker = list(size = 8,
                          opacity = 0.7)) %>%
        layout(title = "Volcano Plot",
               xaxis = list(title = "log2 Fold Change"),
               yaxis = list(title = "-log10 adjusted p-value"))
    })
  })
  
  # Download handler
  output$download_results <- downloadHandler(
    filename = function() {
      input$comparison_select
    },
    content = function(file) {
      write.csv(values$current_results, file, row.names = FALSE)
    }
  )
  
  # Display QC plots
  observe({
    req(values$results_ready)
    
    # RNA degradation plot
    output$rna_deg_plot <- renderPlot({
      plot.new()
      text(0.5, 0.5, "RNA Degradation Plot")
    })
    
    # Box plot
    output$boxplot <- renderPlot({
      plot.new()
      text(0.5, 0.5, "Box Plot")
    })
    
    # MA plot
    output$ma_plot <- renderPlot({
      plot.new()
      text(0.5, 0.5, "MA Plot")
    })
    
    # Density plot
    output$density_plot <- renderPlot({
      plot.new()
      text(0.5, 0.5, "Density Plot")
    })
  })
  
  # Render normalized boxplot
  output$norm_boxplot <- renderPlot({
    req(values$results_ready)
    normalized_data <- read.csv(file.path(input$output_dir, 
                                        "Normalized_Data",
                                        "normalized_expression_values.csv"))
    boxplot(normalized_data,
            main = "Normalized Expression Values",
            las = 2,
            cex.axis = 0.8,
            col = brewer.pal(8, "Set2"))
  })
  
  # Render heatmap
  output$heatmap_plot <- renderPlot({
    req(values$results_ready, input$comparison_select)
    
    # Read the results file
    results_file <- file.path(input$output_dir, "DE_Analysis",
                             input$comparison_select)
    results <- read.csv(results_file)
    
    # Get top 25 genes
    top_genes <- head(order(results$adj.P.Val), 25)
    
    # Get expression data
    expr_data <- as.matrix(read.csv(file.path(input$output_dir,
                                             "Normalized_Data",
                                             "normalized_expression_values.csv"),
                                   row.names = 1))[top_genes, ]
    
    # Get gene symbols
    gene_symbols <- results$Symbol[top_genes]
    
    # Create annotation
    groups <- sapply(seq_along(input$cel_files$name), function(i) {
      input[[paste0("group_", i)]]
    })
    sample_anno <- data.frame(
      Group = groups,
      row.names = colnames(expr_data)
    )
    
    # Create heatmap
    pheatmap(expr_data,
             labels_row = gene_symbols,
             annotation_col = sample_anno,
             scale = "row",
             main = "Top 25 Differentially Expressed Genes",
             fontsize_row = 8,
             fontsize_col = 8)
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server) 