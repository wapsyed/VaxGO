library(shiny)
library(Matrix)
library(circlize)
library(RColorBrewer)
library(org.Hs.eg.db)
library(plotly)
library(clusterProfiler)
library(janitor)
library(ggsci)
library(patchwork)
library(ggthemes)
library(here)
library(readxl)
library(scales)
library(tidyverse)
library(vegan)
library(here)
library(corto)
library(DT)
library(BH)
library(cytolib)
library(heatmaply)
library(BiocManager)
library(devtools)


#Increase upload file size
options(shiny.maxRequestSize=200*1024^2)

################### Custom Functions -------
##### autoORA (for multiple conditions)
autoORA <- function(df,
                    TERM2GENE,
                    pAdjustMethod, 
                    pvalueCutoff_ora, 
                    qvalueCutoff_ora,
                    geneset_name) {
  
  enricher_conditional <- function(condition_2, df, TERM2GENE, pAdjustMethod, pvalueCutoff_ora, qvalueCutoff_ora) {
    genes <- df %>% 
      filter(condition == condition_2) %>% 
      pull(genes)
    
    tryCatch({
      enricher(gene = genes,
               TERM2GENE = TERM2GENE,
               pAdjustMethod = pAdjustMethod,
               pvalueCutoff = pvalueCutoff_ora, 
               qvalueCutoff = qvalueCutoff_ora) %>% 
        as.data.frame() %>%
        mutate(condition = condition_2,
               geneset_name = geneset_name)
    }, error = function(e) {
      NULL
    })
  }
  
  ora_results <- df %>%
    distinct(condition) %>%
    pull(condition) %>%
    map_dfr(~ enricher_conditional(.x, df, TERM2GENE, pAdjustMethod, pvalueCutoff_ora, qvalueCutoff_ora))
  
  ora_results_final <- ora_results %>%
    rownames_to_column("delete") %>%
    dplyr::select(-delete, -Description) %>%
    dplyr::select(condition, geneset_name, everything()) %>% 
    clean_names() %>% 
    rename(process = id)
  
  return(ora_results_final)
}


####### autoGSEA (for multiple conditions)
autoGSEA <- function(df, TERM2GENE, 
                     geneset_name, 
                     minGSSize, 
                     maxGSSize,
                     pvalueCutoff,
                     pAdjustMethod) {
  results <- list()
  
  condicoes <- df %>%
    dplyr::pull(condition) %>%
    unique() %>%
    as.character()
  
  for (condicao in condicoes) {
    degs_condicao <- df %>%
      dplyr::filter(condition == condicao) %>%
      dplyr::select(genes, log2fold_change) %>%
      distinct() %>%
      mutate(rank = rank(log2fold_change, ties.method = "random")) %>%
      arrange(desc(rank))
    
    gene_list_lf2c <- as.vector(degs_condicao$log2fold_change)
    names(gene_list_lf2c) <- degs_condicao$genes
    gene_list_lf2c <- na.omit(gene_list_lf2c)
    gene_list_lf2c <- sort(gene_list_lf2c, decreasing = TRUE)
    
    auto_gsea <- tryCatch({
      GSEA(geneList = gene_list_lf2c,
           TERM2GENE = TERM2GENE,
           minGSSize = minGSSize,
           maxGSSize = maxGSSize,
           pvalueCutoff = pvalueCutoff,
           pAdjustMethod = pAdjustMethod) %>%
        as.data.frame() %>%
        arrange(qvalue) %>%
        mutate(condition = condicao,
               gsea_enrichment = geneset_name)
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(auto_gsea)) {
      results[[paste(condicao, geneset_name, sep = "_")]] <- auto_gsea
    }
  }
  
  final_result <- bind_rows(results) %>%
    rename(process = ID) %>%
    dplyr::select(-Description) %>%
    clean_names() %>%
    dplyr::select(condition, everything())
  
  return(final_result)
}


# Function for overlapping -------
overlap_genes <- function(cond1, cond2, data) {
  genes_cond1 <- data$genes[data$process == cond1]
  genes_cond2 <- data$genes[data$process == cond2]
  
  genes_shared <- intersect(genes_cond1, genes_cond2)
  
  genes_notshared_cond1 <- setdiff(genes_cond1, genes_cond2)
  genes_notshared_cond2 <- setdiff(genes_cond2, genes_cond1)
  
  total_genes_cond1 <- length(genes_cond1)
  total_genes_cond2 <- length(genes_cond2)
  
  percentage_shared_cond1 <- length(genes_shared) / total_genes_cond1 * 100
  percentage_shared_cond2 <- length(genes_shared) / total_genes_cond2 * 100
  
  shared_genes <- data.frame(
    Cond1 = cond1,
    Cond2 = cond2,
    Shared = length(genes_shared),
    NotShared_cond1 = length(genes_notshared_cond1),
    NotShared_cond2 = length(genes_notshared_cond2),
    Total_Genes_Cond1 = total_genes_cond1,
    Total_Genes_Cond2 = total_genes_cond2,
    Genes_Names = paste(genes_shared, collapse = ", "),
    Percentage_Shared_Cond1 = percentage_shared_cond1,
    Percentage_Shared_Cond2 = percentage_shared_cond2
  )
  
  return(shared_genes)
}


# Fisher's exact test
fisher_exact_test <- function(shared, notshared1, notshared2) {
  cont_table <- matrix(c(shared, notshared1, notshared2, 0), nrow = 2)
  results_fisher <- fisher.test(cont_table)
  return(results_fisher$p.value)
}

############################### Load data ------

########## Gene sets ------

######## IMMUNEGO -----
ImmuneGO_Annotated_GO <- readRDS(url("https://github.com/wapsyed/VaxGO/raw/main/Tables/ImmuneGO_Annotated_GO_2024-05-28.rds", method="libcurl"))

ImmuneGO_annotation = ImmuneGO_Annotated_GO %>% 
  dplyr::select(process, immune_system:immune_tissue) 

ImmuneGO_genes_general = ImmuneGO_Annotated_GO %>% 
  separate_rows(genes, sep = ",") %>% 
  filter(go_term == "Manual",
         !process %in% c("ADAPTIVE IMMUNE SYSTEM",
                         "INNATE IMMUNE SYSTEM",
                         "BCR REPERTOIRE",
                         "TCR REPERTOIRE")) %>% 
  dplyr::select(process, genes)

ImmuneGO_genes_specific = ImmuneGO_Annotated_GO %>% 
  separate_rows(genes, sep = ",") %>% 
  filter(!go_term == "Manual") %>% 
  dplyr::select(process=gene_set_short, genes)

######## CELL MARKER IMMUNE ---------

CellMarker_ImmuneCells = readRDS(url("https://github.com/wapsyed/VaxGO/raw/main/Tables/CellMarker_ImmuneCells.rds", method="libcurl")) %>% 
  dplyr::rename(process = cell_name, genes = marker)

CellMarker_annotation = CellMarker_ImmuneCells %>% 
  dplyr::select(process, Type) %>% 
  distinct()

CellMarker_genes = CellMarker_ImmuneCells %>% 
  dplyr::select(-Type)

######## VAX MSIGDB --------
VaxSigDB_Gene_sets_Annotated_RAW <- readRDS(url("https://github.com/wapsyed/VaxGO/raw/main/Tables/VaxSigDB_Gene_sets_Annotated_RAW.rds", method = "libcurl"))

#Filters
vaxsig_columns = c("VACCINE", "PLATFORM", "TARGET_PATHOGEN_DISEASE", "MICROBE_TYPE", "STUDY_TYPE", "STUDY_SUBTYPE","DATE-TIME", "DATE", "TIME", "AGE", "AGE_CATEGORY", "SAMPLE_SOURCE", "SYSTEMATIC_NAME", "GENE_SYMBOLS")
filter_vaxsig_samplesource = "PBMC"
filter_vaxsig_studytype = "VACCINE"
filter_vaxsig_studysubtype = "VAC ONLY"

VaxSigDB_Genesets_filtered = VaxSigDB_Gene_sets_Annotated_RAW %>% 
  dplyr::select(vaxsig_columns) %>% 
  filter(SAMPLE_SOURCE == filter_vaxsig_samplesource,
         STUDY_TYPE == filter_vaxsig_studytype,
         STUDY_SUBTYPE == filter_vaxsig_studysubtype) 
#Annotation
VaxSigDB_annotation = VaxSigDB_Genesets_filtered %>% 
  mutate(process = paste0(VACCINE, " (", `DATE-TIME`, ", ", AGE, " YO)")) %>% 
  dplyr::select(process, PLATFORM, TARGET_PATHOGEN_DISEASE, MICROBE_TYPE, DATE, TIME)
#Genes
VaxSigDB_Genes_filtered = VaxSigDB_Genesets_filtered %>% 
  separate_rows(GENE_SYMBOLS, sep = ",") %>% 
  mutate(process = paste0(VACCINE, " (", `DATE-TIME`, ", ", AGE, " YO)")) %>% 
  dplyr::select(process, genes = GENE_SYMBOLS)

######## BTM MODULES --------
btm_annotation_table <- readRDS(url("https://github.com/wapsyed/VaxGO/raw/main/Tables/btm_annotation_immune_table_labelled.rds", method = "libcurl"))

#Filters
filter_btm_category = "immune"
filter_btm_annotationlevel = c("complete", "partial") #Anyone
# filter_btm_annotationlevel = "complete"
# filter_btm_annotationlevel = "partial"

#Annotation
btm_annotation = btm_annotation_table %>% 
  filter(module_category == filter_btm_category,
         annotation_level %in% filter_btm_annotationlevel) %>% 
  dplyr::select(immune_system, immune_subsystem, cells, composite_name)
#Genes
btm_genes = btm_annotation_table %>% 
  filter(module_category == filter_btm_category,
         annotation_level %in% filter_btm_annotationlevel) %>% 
  dplyr::select(process = composite_name, genes = module_member_genes) %>% 
  separate_rows(genes, sep = ",")

# UI definition
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body {
        font-size: 12px;
      }
      .form-group,
      .shiny-input-container {
        font-size: 12px;
      }
      .shiny-input-container input,
      .shiny-input-container select,
      .shiny-input-container textarea {
        font-size: 12px;
      }
      .btn,
      .btn-primary {
        font-size: 12px;
      }
      .table,
      .table th,
      .table td {
        font-size: 12px;
      }
      /* Ajuste para as abas */
      .nav-tabs > li > a {
        font-size: 12px !important;
      }
      .nav-tabs > li.active > a {
        font-size: 12px !important;
      }
      .nav-tabs > li > a:hover {
        font-size: 12px !important;
      }
      .tab-content {
        font-size: 12px !important;
      }
    "))
  ),
  titlePanel(
    div(
      img(src = "vaxgo_logo_40.png", height = 40, style = "margin-right: 10px;"),
      "VaxGO"
    )
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      textInput(inputId = "filename", 
                label = "Project title", 
                value = "my_condition_degs"),
      
      fileInput("conditions_genes", 
                "Upload your file",
                accept = ".csv"),
      
      h4("Filters", style = "margin-top: 20px;"),
      numericInput(inputId = "filter_padj", 
                   label = "Adjusted p-value cutoff", 
                   value = 0.05, 
                   min = 0, max = 1, step = 0.01),
      
      numericInput(inputId = "filter_logfc", 
                   label = "Log2 Fold Change cutoff", 
                   value = 1, 
                   min = 0, max = 10, step = 0.1),
      actionButton("go_1", "Go"),
      
      h4("GSEA Parameters", style = "margin-top: 20px;"),
      
      numericInput(inputId = "minGSSize", 
                   label = "Minimum Gene Set Size", 
                   value = 1, 
                   min = 1, max = 1000, step = 1),
      
      numericInput(inputId = "maxGSSize", 
                   label = "Maximum Gene Set Size", 
                   value = 1000, 
                   min = 1, max = 10000, step = 1),
      
      selectInput(inputId = "pAdjustMethod", 
                  label = "p-value Adjustment Method", 
                  choices = c("BH", "BY", "fdr", "holm", "hochberg", "hommel", "bonferroni"),
                  selected = "BH"),
      
      numericInput(inputId = "pvalueCutoff_gsea", 
                   label = "GSEA p-value Cutoff", 
                   value = 0.25, 
                   min = 0, max = 1, step = 0.01),
      
      selectInput(inputId = "organism", 
                  label = "Organism", 
                  choices = c("org.Hs.eg.db", "org.Mm.eg.db"),
                  selected = "org.Hs.eg.db"),
      
      actionButton("go_2", "Go"),
      
      h4("ssGSEA Parameters", style = "margin-top: 20px;"),
      fileInput(inputId = "ssgsea_gene", 
                label = "Normalized counts",
                accept = ".csv"),
      fileInput(inputId = "metadata", 
                label = "Sample annotations",
                accept = ".csv"),
      actionButton("go_ssgsea", "Go"),
      selectInput("selected_cellmarker", 
                  "Select a cell:", 
                  choices = NULL,  # Inicialmente vazio, será atualizado
                  selected = NULL,
                  multiple = FALSE),
      
      h4("ORA Parameters", style = "margin-top: 20px;"),
      numericInput(inputId = "minGSSize_ora", 
                   label = "Minimum Gene Set Size", 
                   value = 1, 
                   min = 1, max = 1000, step = 1),
      
      numericInput(inputId = "maxGSSize_ora", 
                   label = "Maximum Gene Set Size", 
                   value = 1000, 
                   min = 1, max = 10000, step = 1),
      
      selectInput(inputId = "pAdjustMethod_ora", 
                  label = "p-value Adjustment Method", 
                  choices = c("BH", "BY", "fdr", "holm", "hochberg", "hommel", "bonferroni"),
                  selected = "BH"),
      
      numericInput(inputId = "pvalueCutoff_ora", 
                   label = "ORA p-value Cutoff", 
                   value = 0.05, 
                   min = 0, max = 1, step = 0.01),
      
      numericInput(inputId = "qvalueCutoff_ora", 
                   label = "ORA q-value Cutoff", 
                   value = 0.10, 
                   min = 0, max = 1, step = 0.01),
      actionButton("go_ora", "Go"),
      
      
      h4("Overlap Parameters", style = "margin-top: 20px;"),
      numericInput(inputId = "filter_padj_overlap", 
                   label = "Adjusted p-value cutoff", 
                   value = 0.05, 
                   min = 0, max = 1, step = 0.01),
      
      numericInput(inputId = "filter_logfc_overlap", 
                   label = "Log2 Fold Change cutoff", 
                   value = 1, 
                   min = 0, max = 10, step = 0.1),
      
      radioButtons(inputId = "overlap_metric", 
                   label = "Overlap method", 
                   choices = c("jaccard_distance", "Shared"),
                   selected = "jaccard_distance"),
      
      actionButton("go_overlap", "Go"),
      
      h4("Correlation Parameters", style = "margin-top: 20px;"),
      numericInput(inputId = "filter_padj_corr", 
                   label = "Adjusted p-value cutoff", 
                   value = 0.05, 
                   min = 0, max = 1, step = 0.01),
      
      radioButtons(inputId = "corr_metric", 
                   label = "Correlation method", 
                   choices = c("Pearson", "Spearman"),
                   selected = "Pearson"),
      
      actionButton("go_corr", "Run Correlation Analysis")
      
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",
        tabPanel("DEGs Analysis",
                 tabsetPanel(
                   type = "tabs",
                   tabPanel("Volcano Plot",
                            plotlyOutput("volcano_plot")),
                   tabPanel("Bar Plot",
                            plotlyOutput("bar_plot"))
                 ),
                 br(),
                 h4("Genes Table"),
                 DT::dataTableOutput("degs_table")
        ),
        tabPanel("GSEA",
                 tabsetPanel(
                   type = "tabs",
                   tabPanel("ImmuneGO General",
                            DT::dataTableOutput("gsea_table_immuneGO_general"),
                            plotlyOutput("gsea_plot_immuneGO_general")),
                   tabPanel("ImmuneGO Specific",
                            DT::dataTableOutput("gsea_table_immuneGO_specific"),
                            plotlyOutput("gsea_plot_immuneGO_specific")),
                   tabPanel("CellMarker Immune",
                            DT::dataTableOutput("gsea_table_cellmarker"),
                            plotlyOutput("gsea_plot_cellmarker")),
                   tabPanel("BTM immune",
                            DT::dataTableOutput("gsea_table_btm_immune"),
                            plotlyOutput("gsea_plot_btm_immune")),
                   tabPanel("Vax MSigDB",
                            DT::dataTableOutput("gsea_table_vaxsigdb"),
                            plotlyOutput("gsea_plot_vaxsigdb"))
                 )
        ),
        tabPanel("ssGSEA",
                 tabsetPanel(
                   type = "tabs",
                   tabPanel("CellMarker",
                            DT::dataTableOutput("ssgsea_table_cellmarker"),
                            plotlyOutput("ssgsea_plot_cellmarker"))
                   
                 )
        ),
        tabPanel("ORA",
                 tabsetPanel(
                   type = "tabs",
                   tabPanel("GeneralGO",
                            DT::dataTableOutput("ora_table_general"),
                            plotOutput("ora_plot_general")),
                   tabPanel("ImmuneGO General",
                            DT::dataTableOutput("ora_table_immuneGO_general"),
                            plotlyOutput("ora_plot_immuneGO_general")),
                   tabPanel("ImmuneGO Specific",
                            DT::dataTableOutput("ora_table_immuneGO_specific"),
                            plotlyOutput("ora_plot_immuneGO_specific")),
                   tabPanel("CellMarker Immune",
                            DT::dataTableOutput("ora_table_cellmarker_immune"),
                            plotlyOutput("ora_plot_cellmarker_immune")),
                   tabPanel("BTM immune",
                            DT::dataTableOutput("ora_table_btm_immune"),
                            plotlyOutput("ora_plot_btm_immune")),
                   tabPanel("VAX MSigDB",
                            DT::dataTableOutput("ora_table_vax_filtered"),
                            plotlyOutput("ora_plot_vax_filtered"))
                 )
        ),
        tabPanel("Gene Overlap",
                 tabsetPanel(
                   type = "tabs",
                   tabPanel("All DEGs",
                            DT::dataTableOutput("overlap_table"),
                            plotlyOutput("overlap_heatmap")),
                   tabPanel("Immune",
                            DT::dataTableOutput("overlap_table_immune"),
                            plotlyOutput("overlap_heatmap_immune")),
                   tabPanel("Not immune",
                            DT::dataTableOutput("overlap_table_not_immune"),
                            plotlyOutput("overlap_heatmap_not_immune"))
                   
                 )
        ),
        tabPanel("Correlation",
                 tabsetPanel(
                   type = "tabs",
                   tabPanel("All DEGs",
                            DT::dataTableOutput("corr_table")
                   )
                   # tabPanel("Immune",
                   #          DT::dataTableOutput("overlap_table_immune"),
                   #          plotlyOutput("overlap_heatmap_immune")),
                   # tabPanel("Not immune",
                   #          DT::dataTableOutput("overlap_table_not_immune"),
                   #          plotlyOutput("overlap_heatmap_not_immune"))
                   
                 )
        )
      )
    )
  )
)

# Server definition
server <- function(input, output, session) {
  
  # Reactive expression to load and process the file
  conditions_genes <- reactive({
    req(input$conditions_genes)
    
    # Load the file
    conditions_genes <- read.csv(input$conditions_genes$datapath) %>% 
      distinct()
    conditions_genes
  })
  
  # Reactive expression to filter significant DEGs
  degs_sig <- reactive({
    req(conditions_genes(), input$filter_padj)
    
    degs_sig <- conditions_genes() %>%
      filter(padj < input$filter_padj) %>% 
      distinct()
    
    degs_sig
  })
  
  # Reactive expression to further filter DEGs and extract gene names
  genes <- reactive({
    req(degs_sig())
    
    genes <- degs_sig() %>%
      filter(log2fold_change >= input$filter_logfc | log2fold_change <= - input$filter_logfc,
             padj < input$filter_padj) %>%
      pull(genes)
    
    genes
  })
  
  #Volcano plot
  output$volcano_plot <- renderPlotly({
    req(input$go_1) # Action button
    req(conditions_genes())
    
    # Prepare the data
    data <- conditions_genes() %>%
      mutate(log_q = -log10(padj))
    
    # List of unique conditions
    unique_conditions <- unique(data$condition)
    
    # Create an empty plot
    p <- plot_ly()
    
    # Add a trace for each condition
    for (cond in unique_conditions) {
      p <- p %>%
        add_trace(
          data = data[data$condition == cond, ],
          x = ~log2fold_change, 
          y = ~log_q, 
          text = ~genes,
          type = 'scatter', 
          mode = 'markers',
          marker = list(
            color = ~log2fold_change,
            size = 10,
            opacity = 0.9,
            colorscale = c("#3399FF", "white", "#FF6666")
          ),
          name = cond,
          visible = ifelse(cond == unique_conditions[1], TRUE, FALSE) # Show only the first condition by default
        )
    }
    
    # Create updatemenus based on conditions
    condition_buttons <- lapply(1:length(unique_conditions), function(i) {
      visibility <- rep(FALSE, length(unique_conditions))
      visibility[i] <- TRUE
      list(
        method = "restyle",
        args = list("visible", visibility),
        label = unique_conditions[i]
      )
    })
    
    # Add layout options including the update menus
    p %>%
      layout(
        title = "Volcano Plot",
        xaxis = list(title = "Log2 Fold Change"),
        yaxis = list(title = "-Log10 Adjusted P-value"),
        updatemenus = list(
          list(
            y = 1,  # Position the buttons near the top
            x = 0,   # Center the buttons horizontally
            xanchor = "center",
            yanchor = "top",
            buttons = condition_buttons,
            showactive = TRUE
          )
        )
      )
  })
  
  # Bar Plot
  output$bar_plot <- renderPlotly({
    req(input$go_1) # Action button
    req(conditions_genes())
    
    # Prepare the data
    data <- conditions_genes()  %>%
      filter(log2fold_change >= input$filter_logfc | log2fold_change <= - input$filter_logfc,
             padj < input$filter_padj) %>% 
      group_by(condition, direction) %>% 
      summarize(n_direct = n()) %>% 
      ungroup()
    
    
    # Create the volcano plot directly with plotly
    plot_ly(data, 
            x = ~condition, 
            y = ~n_direct, 
            color = ~direction, 
            colors = c("#3399FF", "#FF6666"),
            text = ~direction,
            type = 'bar') %>%
      layout(
        title = "DEGs Counts",
        xaxis = list(title = "Condition"),
        yaxis = list(title = "# Genes"),
        barmode = 'group',  
        coloraxis = list(colorbar = list(title = 'Direction'))
      )
  })
  
  #DEGs table
  output$degs_table <- DT::renderDataTable({
    req(input$go_1) # Action button
    req(degs_sig())
    
    # Define input values for the cutoffs
    filter_l2fc <- input$filter_logfc
    filter_neg_l2fc <- -input$filter_logfc
    
    # Transform the filtered data into a DataTable
    degs_sig() %>%
      datatable(options = list(pageLength = 25)) %>%
      formatStyle(
        'log2fold_change',
        backgroundColor = styleInterval(
          c(filter_neg_l2fc, filter_l2fc),
          c('#3399FF', 'white', '#FF6666')
        ),
        fontWeight = 'bold'  # Optional: make the text bold for better visibility
      )
  })
  
  #GSEA ImmuneGO general
  # Reactive expression to perform GSEA analyses
  gsea_results_ImmuneGO_General <- reactive({
    req(input$go_2)
    req(degs_sig())
    
    degs_sig = degs_sig()
    
    # Perform ORA for ImmuneGO general
    degs_GSEA_ImmuneGO_General <- degs_sig %>%
      dplyr::select(condition, genes, log2fold_change) %>% 
      autoGSEA(TERM2GENE = ImmuneGO_genes_general, 
               geneset_name = "ImmuneGO General",
               minGSSize = input$minGSSize,
               maxGSSize = input$maxGSSize,
               pvalueCutoff = input$pvalueCutoff_gsea,
               pAdjustMethod = input$pAdjustMethod) 
  })
  
  # Render the DataTable
  output$gsea_table_immuneGO_general <- DT::renderDataTable({
    req(input$go_2)
    req(gsea_results_ImmuneGO_General())
    gsea_results_ImmuneGO_General() %>% 
      datatable(extensions = 'Buttons',
                
                options = list(
                  paging = TRUE,
                  searching = TRUE,
                  fixedColumns = TRUE,
                  autoWidth = TRUE,
                  ordering = TRUE,
                  dom = 'Blfrtip',
                  buttons = c('copy', 'csv', 'excel')
                ),
                
                class = "display"
      )
  })
  #GSEA ImmuneGO general Plot
  output$gsea_plot_immuneGO_general <- renderPlotly({
    req(input$go_2) # Action button
    req(gsea_results_ImmuneGO_General())
    
    # Prepare the data
    data <- gsea_results_ImmuneGO_General() %>% 
      group_by(condition) %>%
      arrange(condition, desc(-log(qvalue))) %>% 
      slice_head(n = 5) %>%
      ungroup() %>%
      mutate(process = fct_reorder(process, -log(qvalue)))  # Reordenar o processo
    
    num_conditions <- data %>%
      distinct(condition) %>%
      nrow()
    
    ggplot(data) +
      aes(x = -log(qvalue), y = process,
          fill = -log(qvalue)) +
      geom_col() +
      facet_wrap(~condition, scales = "free_y", nrow = num_conditions) +
      labs(x = "-log10(q-value)", y = "Process", title = "GSEA Results") +
      theme_minimal()
  })
  
  
  
  #GSEA ImmuneGO specific
  # Reactive expression to perform GSEA analyses
  gsea_results_ImmuneGO_genes_specific <- reactive({
    req(input$go_2)
    req(degs_sig())
    
    degs_sig = degs_sig()
    
    # Perform ORA for ImmuneGO general
    degs_GSEA_ImmuneGO_genes_specific <- degs_sig %>%
      dplyr::select(condition, genes, log2fold_change) %>% 
      autoGSEA(TERM2GENE = ImmuneGO_genes_specific, 
               geneset_name = "ImmuneGO Specific",
               minGSSize = input$minGSSize,
               maxGSSize = input$maxGSSize,
               pvalueCutoff = input$pvalueCutoff_gsea,
               pAdjustMethod = input$pAdjustMethod) 
  })
  
  # Render the DataTable
  output$gsea_table_immuneGO_specific <- DT::renderDataTable({
    req(input$go_2)
    req(gsea_results_ImmuneGO_genes_specific())
    gsea_results_ImmuneGO_genes_specific() %>% 
      datatable(extensions = 'Buttons',
                
                options = list(
                  paging = TRUE,
                  searching = TRUE,
                  fixedColumns = TRUE,
                  autoWidth = TRUE,
                  ordering = TRUE,
                  dom = 'Blfrtip',
                  buttons = c('copy', 'csv', 'excel')
                ),
                
                class = "display"
      )
  })
  #GSEA ImmuneGO specific Plot
  output$gsea_plot_immuneGO_specific <- renderPlotly({
    req(input$go_2) # Action button
    req(gsea_results_ImmuneGO_genes_specific())
    
    # Prepare the data
    data <- gsea_results_ImmuneGO_genes_specific()%>% 
      group_by(condition) %>%
      arrange(condition, desc(-log(qvalue))) %>% 
      slice_head(n = 5) %>%
      ungroup() %>%
      mutate(process = fct_reorder(process, -log(qvalue)))  # Reordenar o processo
    
    num_conditions <- data %>%
      distinct(condition) %>%
      nrow()
    
    ggplot(data) +
      aes(x = -log(qvalue), y = process,
          fill = -log(qvalue)) +
      geom_col() +
      facet_wrap(~condition, scales = "free_y", nrow = num_conditions) +
      labs(x = "-log10(q-value)", y = "Process", title = "GSEA Results") +
      theme_minimal()})