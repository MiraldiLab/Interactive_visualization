library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(circlize)


counts <- read.delim('/Volumes/miraldiNB/anthony/projects/HAE/analysis/230209_pseudobulk_scrna_major/scrna_IFN_bulk_major_group_combat.txt',sep='\t',header=T,row.names = 1)


ui = fluidPage(
  h3("My first interactive ComplexHeatmap Shiny app"),
  p("This is an interactive heatmap visualization on a random matrix."),
  # Sidebar layout with a input and output definitions
  sidebarLayout (
    # Inputs: Select variables to plot
    sidebarPanel (
      #Subset of genes to plot
      fileInput("genes_to_plot", "Select a .txt file of genes to plot (optional)"),
      
      # Select variable for Row Normalization
      selectInput (inputId = "norm_row", label = "Normalization (rows)",
                   choices = c ("z-score_across", "z-score_within", "log2fc_mean"),
                   selected = "z-score_across"),
      
      # Select variable for Number of Clusters
      sliderInput(inputId = "row_kmn", label = "Row Clusters", 
                  value = 1, min = 1, max = 10),
      
      #color range slider
      sliderInput("colbar_range", label = h3("Color Bar Range"), min = -10, 
                  max = 10, value = c(-2, 2)),
      
      #Button to generate heatmap
      actionButton("generate_heatmap", "Generate Heatmap")
    ),
    # Output: Show scatterplot
    mainPanel (
      InteractiveComplexHeatmapOutput(),
    )
  )
)

server = function(input, output, session) {
  observeEvent(input$generate_heatmap,{
    if (is.null(input$genes_to_plot)) {
      rows <- rownames(counts[1:1000,])
    } else {
      rows <- readLines(input$genes_to_plot$datapath)
    }
    
    counts <- counts[rows,]
    
    if(input$norm_row == 'z-score_across'){
      counts_norm <- t(scale(t(counts)))
    } else if (input$norm_row == 'log2fc_mean'){
      rowmean <- rowMeans(counts)
      counts$rowmean <- rowmean
      counts_norm <- log2(counts/counts$rowmean)
      counts_norm <- counts_norm[,-ncol(counts_norm)]
    } else {
      counts_norm <- counts
    }
    
    heat_col <- colorRamp2(c(input$colbar_range[1], 0, input$colbar_range[2]), c('dodgerblue2','white','red'))

    ht = Heatmap(counts_norm[rows,sort(colnames(counts_norm))],
                 row_km = as.numeric(input$row_kmn), 
                 cluster_columns = FALSE,
                 column_names_side = 'top',
                 row_names_side = 'left',
                 col=heat_col)
    makeInteractiveComplexHeatmap(input, output, session, ht)
  })
}

shinyApp(ui, server)
