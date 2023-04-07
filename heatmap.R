library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(circlize)
library(igraph)
library(shiny)
library(visNetwork)

ui <- fluidPage(
  h3("Modeling the heterogeneous interferon response of the human airway epithelium at single cell resolution"),
  p("Citation..."),
  
  #Tab 1: Gene Expression Visualization
  tabsetPanel(
    tabPanel("Gene Expression", fluid=T,
      p("Interactive heatmap visualization of the pseudobulk gene expression matrix."),
      sidebarLayout (
        sidebarPanel (
          #Gene Expression
          fileInput("gene_expression", "Select a Gene Expression file (Required)"),
          
          #Metadata
          fileInput("gene_expression_metadata", "Select a Metadata file (Optional)"),
          
          #Subset of genes to plot
          fileInput("genes_to_plot", "Select a .txt file of genes to plot (optional)"),
          
          # Select variable for Row Normalization
          selectInput (inputId = "norm_row", label = "Normalization (rows)",
                       choices = c ("None","z-score across", "log2fc (mean)"),
                       selected = "none"),
          
          # Select variable for Number of Clusters
          sliderInput(inputId = "row_kmn", label = "Row Clusters", 
                      value = 1, min = 1, max = 10),
          
          #color range slider
          sliderInput("colbar_range", label = "Color Bar Range", min = -10, 
                      max = 10, value = c(-2, 2)),
          
          #Button to generate heatmap
          actionButton("generate_heatmap_gene_exp", "Generate Heatmap")
        ),
        
        # Output: Show Heatmap
        mainPanel (
          InteractiveComplexHeatmapOutput(heatmap_id = 'ht_gene_exp', height1 = 700,width1=500,height2 = 700,width2=500),
        )
      )
    ),
    
  #Tab 2: GSEA Results Visualization
  tabPanel("GSEA (Under construction)", fluid=T,
            p("Interactive visualization of GSEA results."),      
               # Sidebar layout with a input and output definitions
               sidebarLayout (
                 sidebarPanel (
                   #GSEA Results
                   fileInput("gsea_res", "Select a GSEA result file"),
                   
                   #padj cutodd
                   textInput(inputId='padj_cutoff_gsea', label='FDR cutoff', value = "0.1", placeholder = 'Default = 0.1'),
                   
                   #Button to generate heatmap
                   actionButton("generate_heatmap", "Generate Heatmap")
                  ),
                 
                # Output: Show Heatmap
                mainPanel (
                  InteractiveComplexHeatmapOutput(),
                )
              )
            ),
  
  #Tab 3: Network Visualization
  tabPanel("Network", fluid=T,
           p("Interactive visualization of Gene Regulatory Networks."),
           sidebarLayout (
             # Inputs: Select variables to plot
             sidebarPanel (
               #Subset of genes to plot
               fileInput("network_to_plot", "Select a network file to plot"),
               
               #Gene Expression
               fileInput("gene_expression_network", "Select a Gene Expression file"),
               
               # Select variable for Row Normalization
               selectInput(inputId = "gene_expression_color", label = "Condition to Color by",
                           choices = NULL,
                           selected = NULL),
               # Select variable for Row Normalization of gene expression
               selectInput (inputId = "norm_row_network", label = "Normalization of Gene Expression",
                            choices = c ("none","z-score", "log2fc"),
                            selected = "none"),
               
               #color range slider
               sliderInput("colbar_range", label = "Node Color Range", min = -10, 
                           max = 10, value = c(-2, 2)),
               
               # Select variable cutoff threshold
               sliderInput(inputId = "threshold_cutoff", label = "Threshold (Quantile)", 
                           value = 0, min = 0, max = 1,step=0.1),
               
               #TFs to draw
               textInput(inputId='tfs_to_plot', label='TFs to plot (space separated)', value = "", placeholder = 'Enter TFs'),
               
               #TF text size
               textInput(inputId='tf_text_size', label='TF text size', value = 50, placeholder = 'Default = 50'),
               
               #Target text size
               textInput(inputId='target_text_size', label='Target text size', value = 50, placeholder = 'Default = 50'),
               
               #Button to generate heatmap
               actionButton("draw_network", "Draw Network")
             ),
             
             mainPanel(
               visNetworkOutput("distPlot", width='100%'),
               textOutput(outputId = "selected_var")
             )
           )
          )
  )
)

server <- function(input, output, session) {
  options(shiny.maxRequestSize=100*1024^2)
  
  #Tab 1: Gene Expression Output
  observeEvent(input$generate_heatmap_gene_exp,{
    
    counts <- read.delim(input$gene_expression$datapath, sep='\t',header=T,row.names=1)
    
    # if (is.null(input$gene_expression_metadata)) {
    #   metadata <- read.delim(input$gene_expression_metadata$datapath, sep='\t',header=T,row.names=1)
    # }
    
    if (is.null(input$genes_to_plot)) {
      rows <- rownames(counts[1:1000,])
    } else {
      rows <- readLines(input$genes_to_plot$datapath)
      #rows <- rows[which(rows %in% rownames(counts))]
    }
    
    counts <- counts[rows,]
    
    if(input$norm_row == 'z-score across'){
      counts_norm <- t(scale(t(counts)))
    } else if (input$norm_row == 'log2fc mean'){
      rowmean <- rowMeans(counts)
      counts$rowmean <- rowmean
      counts_norm <- log2(counts/counts$rowmean)
      counts_norm <- counts_norm[,-ncol(counts_norm)]
    } else {
      counts_norm <- counts
    }
    
    heat_col <- colorRamp2(c(input$colbar_range[1], 0, input$colbar_range[2]), c('dodgerblue2','white','red'))

    ht_gene_exp = Heatmap(counts_norm[rows,sort(colnames(counts_norm))],
                 row_km = as.numeric(input$row_kmn), 
                 cluster_columns = FALSE,
                 column_names_side = 'top',
                 row_names_side = 'left',
                 col=heat_col,
                 show_row_names=F)
    ht_gene_exp = draw(ht_gene_exp)
      makeInteractiveComplexHeatmap(input, output, session, ht_gene_exp, heatmap_id = 'ht_gene_exp')
  })
  
  #Tab 2: GSEA Output
  observeEvent(input$generate_heatmap_gsea,{
    counts_gsea <- read.delim(input$gsea_res$datapath, sep='\t',header=T)
    
  })
  
  #Tab 3: Network Output
  values <- reactiveValues(df_data = NULL)
  y_vals <- NULL
  
  observeEvent(input$gene_expression_network, {
    values$df_data <- read.delim(input$gene_expression_network$datapath, sep='\t',header=T,row.names=1)
    y_vals <- sort(colnames(values$df_data))
    updateSelectInput(
      session = session, 
      inputId = "gene_expression_color",
      choices = y_vals,
      selected = head(y_vals, 1)
    )
    
  })
  
  observeEvent(input$draw_network,{
    if (input$norm_row_network == 'none'){
      values$df_data <- values$df_data
    } else if (input$norm_row_network == 'z-score'){
      values$df_data <- as.data.frame(t(scale(t(values$df_data))))
    } else {
      rowmean <- rowMeans(values$df_data)
      values$df_data$rowmean <- rowmean
      values$df_data <- log2(values$df_data/values$df_data$rowmean)
      values$df_data <- values$df_data[,-ncol(values$df_data)]
    }
    gene_expression_cols <- colorRamp2(c(input$colbar_range[1], 0, input$colbar_range[2]), c('blue','white','red'))
    
    network <- read.delim(input$network_to_plot$datapath, sep='\t', header=T) #read network
    tfs <- unique(network$TF)
    
    if (length(input$tfs_to_plot) == 0) {
      tfs_to_plot <- tfs
      output$selected_var <- renderText({ 
        'Plotting all Tfs'})
    } else if (length(scan(text = input$tfs_to_plot, what = "")) == 1 & input$tfs_to_plot %in% tfs) {
      tfs_to_plot <- input$tfs_to_plot
    } else if (length(scan(text = input$tfs_to_plot, what = "")) == 1 & length(setdiff(input$tfs_to_plot, tfs)) > 0) {
      tfs_to_plot <- tfs
      output$selected_var <- renderText({ 
        paste(input$tfs_to_plot, 'not in network. Defaulting to all TFs.',sep=' ')})
    } else {
      tfs_to_plot <- scan(text = input$tfs_to_plot, what = "")
    }
    
    if(length(setdiff(input$tfs_to_plot,network$TF)) > 1){
      output$selected_var <- renderText({ 
        paste(setdiff(tfs_to_plot,network$TF), 'not in network. Removing from input.',sep=' ')})
      tfs_to_plot <- setdiff(tfs_to_plot, setdiff(input$tfs_to_plot,network$TF))
    }
    
    network <- network[which(network$TF %in% tfs_to_plot),] #subset by TF
    network <- subset(network, combStability >= quantile(combStability, input$threshold_cutoff))
    
    tfs <- unique(network$TF)
    targets <- unique(network$Target)
    nodes <- data.frame(id=union(tfs, targets))
    nodes$label <- nodes$id
    nodes$shape <- ifelse(nodes$label %in% tfs, 'diamond','dot')
    nodes$physics <- FALSE
    nodes$size <- 50
    nodes$font.size <- ifelse(nodes$label %in% tfs, input$tf_text_size, input$target_text_size)
    
    node_colors <- values$df_data[nodes$id,input$gene_expression_color]
    node_colors <- gene_expression_cols(node_colors)
    
    nodes$color <- gsub('{2}$','',node_colors)
    
    edges <- data.frame(from=network$TF,
                        to=network$Target,
                        arrows='to',
                        title='',
                        value=abs(network$SignedQuantile),
                        color=network$stroke,
                        dashes=ifelse(network$stroke.dasharray =='None', F, T),
                        physics=FALSE)
    
    output$distPlot <- renderVisNetwork({
      set.seed(125)
      
      visNetwork(nodes, edges,width='100%', height='100%') %>%
        visOptions(highlightNearest = list(enabled =FALSE, degree = 1, hideColor='rgba(255,255,255,0.5)'))%>%
        visEdges(smooth = FALSE,arrows = "middle")%>%
        visPhysics(enabled = FALSE, stabilization = FALSE)%>%
        visIgraphLayout(type = "full")
      # %>% 
      #   visExport(type = "png", name = "export-network", 
      #             float = "left", label = "Save network", style= "") 
      
    })
  })
}

shinyApp(ui, server)
