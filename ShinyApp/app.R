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
      print(paste("Error in GSEA for condition:", condicao, "Error:", e$message))
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
        background-color: #ffffff;
        color: #000000;
        font-family: 'Montserrat', sans-serif;
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
      pre, code {
        font-family: 'Montserrat', monospace;
      }
    ")),
    tags$link(
      href = "https://fonts.googleapis.com/css2?family=Montserrat:wght@400;700&display=swap",
      rel = "stylesheet"
    ),
    tags$link(
      rel = "stylesheet",
      href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css"
    ),
    tags$link(
      rel = "stylesheet",
      href = "https://code.ionicframework.com/ionicons/2.0.1/css/ionicons.min.css"
    )
  ),
  div(
    style = "display: flex; justify-content: space-between; align-items: center;",
    titlePanel(
      div(
        h2("VaxGO Tool", style = "margin-bottom: 0px;"),
        tags$p(
          HTML(
            paste0(
              "Developer: Wasim Aluísio Prates-Syed ",
              "<a href='https://github.com/wapsyed' target='_blank'>",
              "github.com/wapsyed",
              "</a>"
            )
          ),
          style = "font-size: 14px; margin-top: 5px;"
        )
      )
    ),
    tags$nav(
      style = "display: flex; align-items: center;",
      tags$a(
        href = "https://instagram.com/wasimvacinas",
        target = "_blank",
        style = "margin-left: 10px;",
        tags$i(class = "fa fa-instagram", style = "font-size: 20px;")
      ),
      tags$a(
        href = "https://www.linkedin.com/in/wapsyed",
        target = "_blank",
        style = "margin-left: 10px;",
        tags$i(class = "fa fa-linkedin", style = "font-size: 20px;")
      ),
      tags$a(
        href = "https://github.com/wapsyed",
        target = "_blank",
        style = "margin-left: 10px;",
        tags$i(class = "fa fa-github", style = "font-size: 20px;")
      )
    )
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Parameters"), # Title for the sidebar
      tabsetPanel(
        type = "tabs",
        tabPanel("DEGs",
                 tags$p(
                   HTML(
                     paste0(
                       "Download an ",
                       "<a href='https://github.com/wapsyed/VaxGO/blob/662205136a565e3420c26bfce665073641c45cc4/Example/Example_DEGs_BNT_multiple.csv'>",
                       "example",
                       "</a>"
                     )
                   ),
                 ),
                 
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
                 actionButton("go_1", "Go")
        ),
        tabPanel("GSEA",
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
                 
                 actionButton("go_2", "Go")
        ),
        tabPanel("ssGSEA",
                 h4("ssGSEA Parameters", style = "margin-top: 20px;"),
                 
                 
                 tags$p(
                   HTML(
                     paste0(
                       "Download an ",
                       "<a href='https://github.com/wapsyed/VaxGO/blob/662205136a565e3420c26bfce665073641c45cc4/Example/Example_normalized_counts.csv'>",
                       "example"
                     )
                   ),
                 ),
                 fileInput(inputId = "ssgsea_gene", 
                           label = "Normalized counts",
                           accept = ".csv"),
                 
                 
                 tags$p(
                   HTML(
                     paste0(
                       "Download an ",
                       "<a href='https://github.com/wapsyed/VaxGO/blob/662205136a565e3420c26bfce665073641c45cc4/Example/Example_annotations_samples.csv'>",
                       "example",
                       "</a>"
                     )
                   ),
                 ),
                 fileInput(inputId = "metadata", 
                           label = "Sample annotations",
                           accept = ".csv"),
                 
                 actionButton("go_ssgsea", "Go"),
                 selectInput("selected_cellmarker", 
                             "Select a cell:", 
                             choices = NULL,  # Inicialmente vazio, será atualizado
                             selected = NULL,
                             multiple = FALSE)
        ),
        tabPanel("Overlap",
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
                 
                 actionButton("go_overlap", "Go")
        )
      )
    ),
    mainPanel(
      width = 9,
      h4("Results"), # Title for the sidebar
      
      tabsetPanel(
        type = "tabs",
        tabPanel("DEGs Analysis",
                 fluidRow(
                   column(width = 6,
                          plotlyOutput("volcano_plot")
                   ),
                   column(width = 6,
                          plotlyOutput("bar_plot")
                   )
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
                 DT::dataTableOutput("ssgsea_table_cellmarker"),
                 plotlyOutput("ssgsea_plot_cellmarker")
        ),
        tabPanel("Gene Overlap",
                 plotlyOutput("overlap_heatmap"),
                 DT::dataTableOutput("overlap_table"),
                 br(),
                 plotlyOutput("overlap_heatmap_immune"),
                 DT::dataTableOutput("overlap_table_immune"),
                 br(),
                 plotlyOutput("overlap_heatmap_not_immune"),
                 DT::dataTableOutput("overlap_table_not_immune"),
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
    
    degs_sig_data <- degs_sig()
    print("Data for GSEA ImmuneGO General:")
    print(head(degs_sig_data))
    
    # Perform ORA for ImmuneGO general
    degs_GSEA_ImmuneGO_General <- degs_sig_data %>%
      dplyr::select(condition, genes, log2fold_change) %>% 
      autoGSEA(TERM2GENE = ImmuneGO_genes_general, 
               geneset_name = "ImmuneGO General",
               minGSSize = input$minGSSize,
               maxGSSize = input$maxGSSize,
               pvalueCutoff = input$pvalueCutoff_gsea,
               pAdjustMethod = input$pAdjustMethod) 
    
    print("GSEA ImmuneGO General results:")
    print(head(degs_GSEA_ImmuneGO_General))
    
    degs_GSEA_ImmuneGO_General
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
    
    degs_sig_data <- degs_sig()
    print("Data for GSEA ImmuneGO Specific:")
    print(head(degs_sig_data))
    
    # Perform ORA for ImmuneGO general
    degs_GSEA_ImmuneGO_genes_specific <- degs_sig_data %>%
      dplyr::select(condition, genes, log2fold_change) %>% 
      autoGSEA(TERM2GENE = ImmuneGO_genes_specific, 
               geneset_name = "ImmuneGO Specific",
               minGSSize = input$minGSSize,
               maxGSSize = input$maxGSSize,
               pvalueCutoff = input$pvalueCutoff_gsea,
               pAdjustMethod = input$pAdjustMethod) 
    
    print("GSEA ImmuneGO Specific results:")
    print(head(degs_GSEA_ImmuneGO_genes_specific))
    
    degs_GSEA_ImmuneGO_genes_specific
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
      theme_minimal()
  })
  
  #GSEA CellMarker
  # Reactive expression to perform GSEA analyses
  gsea_results_CellMarker <- reactive({
    req(input$go_2)
    req(degs_sig())
    
    degs_sig_data <- degs_sig()
    print("Data for GSEA CellMarker:")
    print(head(degs_sig_data))
    
    # Perform ORA for ImmuneGO general
    degs_GSEA_CellMarker <- degs_sig_data %>%
      dplyr::select(condition, genes, log2fold_change) %>% 
      autoGSEA(TERM2GENE = CellMarker_genes, 
               geneset_name = "CellMarker immune",
               minGSSize = input$minGSSize,
               maxGSSize = input$maxGSSize,
               pvalueCutoff = input$pvalueCutoff_gsea,
               pAdjustMethod = input$pAdjustMethod) 
    
    print("GSEA CellMarker results:")
    print(head(degs_GSEA_CellMarker))
    
    degs_GSEA_CellMarker
  })
  
  # Render the DataTable
  output$gsea_table_cellmarker <- DT::renderDataTable({
    req(input$go_2)
    req(gsea_results_CellMarker())
    gsea_results_CellMarker() %>% 
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
  #GSEA CellMarker Plot
  output$gsea_plot_cellmarker <- renderPlotly({
    req(input$go_2) # Action button
    req(gsea_results_CellMarker())
    
    # Prepare the data
    data <- gsea_results_CellMarker()%>% 
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
  
  
  #GSEA BTM immune
  # Reactive expression to perform GSEA analyses
  gsea_results_BTM_Immune  <- reactive({
    req(input$go_2)
    req(degs_sig()) 
    
    degs_sig_data <- degs_sig()
    print("Data for GSEA BTM Immune:")
    print(head(degs_sig_data))
    
    # Perform ORA for ImmuneGO general
    degs_GSEA_BTM_Immune <- degs_sig_data %>%
      dplyr::select(condition, genes, log2fold_change) %>% 
      autoGSEA(TERM2GENE = btm_genes, 
               geneset_name = "BTM Immune",
               minGSSize = input$minGSSize,
               maxGSSize = input$maxGSSize,
               pvalueCutoff = input$pvalueCutoff_gsea,
               pAdjustMethod = input$pAdjustMethod) 
    
    print("GSEA BTM Immune results:")
    print(head(degs_GSEA_BTM_Immune))
    
    degs_GSEA_BTM_Immune
  })
  
  # Render the DataTable
  output$gsea_table_btm_immune <- DT::renderDataTable({
    req(input$go_2)
    req(gsea_results_BTM_Immune())
    gsea_results_BTM_Immune() %>% 
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
  #GSEA BTM immune Plot
  output$gsea_plot_btm_immune <- renderPlotly({
    req(input$go_2) # Action button
    req(gsea_results_BTM_Immune())
    
    # Prepare the data
    data <- gsea_results_BTM_Immune()%>% 
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
  
  
  #GSEA VaxSigDB
  # Reactive expression to perform GSEA analyses
  gsea_results_VaxSigDB_Genes_filtered <- reactive({
    req(input$go_2)
    req(degs_sig())
    
    degs_sig_data <- degs_sig()
    print("Data for GSEA VaxSigDB:")
    print(head(degs_sig_data))
    
    # Perform ORA for ImmuneGO general
    degs_GSEA_VaxSigDB_Genes_filtered <- degs_sig_data %>%
      dplyr::select(condition, genes, log2fold_change) %>% 
      autoGSEA(TERM2GENE = VaxSigDB_Genes_filtered, 
               geneset_name = "VaxSigDB Filtered",
               minGSSize = input$minGSSize,
               maxGSSize = input$maxGSSize,
               pvalueCutoff = input$pvalueCutoff_gsea,
               pAdjustMethod = input$pAdjustMethod) 
    
    print("GSEA VaxSigDB results:")
    print(head(degs_GSEA_VaxSigDB_Genes_filtered))
    
    degs_GSEA_VaxSigDB_Genes_filtered
  })
  
  # Render the DataTable
  output$gsea_table_vaxsigdb <- DT::renderDataTable({
    req(input$go_2)
    req(gsea_results_VaxSigDB_Genes_filtered())
    gsea_results_VaxSigDB_Genes_filtered() %>% 
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
  #GSEA VaxSigDB Plot
  output$gsea_plot_vaxsigdb <- renderPlotly({
    req(input$go_2) # Action button
    req(gsea_results_VaxSigDB_Genes_filtered())
    
    # Prepare the data
    data <- gsea_results_VaxSigDB_Genes_filtered()%>% 
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
  
  
  
  #ssGSEA CellMarker
  # Carregamento e processamento dos dados
  ssgsea_results_CellMarker <- reactive({
    req(input$go_ssgsea)  
    req(input$ssgsea_gene, input$metadata)  
    
    # Ler os arquivos carregados
    expression_data <- read.csv(input$ssgsea_gene$datapath) %>% 
      column_to_rownames("genes")
    
    metadata <- read.csv(input$metadata$datapath)
    
    # Obter os nomes das amostras
    sample_names <- colnames(expression_data) %>% 
      as.data.frame() %>% 
      rename(sample_names = '.') %>% 
      distinct()
    
    # ImmuneGO General
    genelist_cellmarker <- CellMarker_ImmuneCells %>%
      dplyr::select(process, genes)
    
    genelist_cellmarker = split(genelist_cellmarker$genes,
                                genelist_cellmarker$process)
    
    
    # Calcular ssGSEA
    nesmat_result_cellmarker <- ssgsea(inmat = expression_data, 
                                       groups = genelist_cellmarker) %>% 
      as.data.frame() %>% 
      t() %>% 
      as.data.frame() %>% 
      rownames_to_column("sample") %>% 
      relocate(sample, .before = everything()) %>% 
      cbind(sample_names) %>% 
      dplyr::select(-sample) %>% 
      column_to_rownames("sample_names") %>% 
      t() %>% 
      as.data.frame() %>% 
      rownames_to_column("process") %>% 
      relocate(process, .before = everything()) %>% 
      pivot_longer(cols = -process, names_to = "sample", values_to = "nes") %>%
      mutate(pvalue = z2p(nes), # Converter NES para p-value
             qvalue = p.adjust(pvalue, method = "BH"),
             logq = -log10(qvalue),
             geneset_name = "CellMarker Immune") %>%
      inner_join(metadata, by = "sample") %>%
      inner_join(CellMarker_ImmuneCells, by = "process") %>%
      dplyr::select(-genes) %>%
      distinct() 
    
    nesmat_result_cellmarker %>% as.data.frame()
  })
  
  # Renderizar o DataTable
  output$ssgsea_table_cellmarker <- DT::renderDataTable({
    req(input$go_ssgsea)
    req(ssgsea_results_CellMarker())
    
    ssgsea_results_CellMarker() %>% 
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
  
  
  observe({
    results <- ssgsea_results_CellMarker()
    if (!is.null(results) && nrow(results) > 0) {
      updateSelectInput(session, "selected_cellmarker", 
                        choices = unique(results$process))
    }
  })
  
  #ssGSEA cellmarker Plot
  output$ssgsea_plot_cellmarker <- renderPlotly({
    req(input$go_ssgsea) 
    req(input$selected_cellmarker)
    req(ssgsea_results_CellMarker())
    
    # Prepare the data
    data <- ssgsea_results_CellMarker() %>% 
      filter(process == input$selected_cellmarker) %>% 
      ggplot(.) +
      aes(x = condition, y = nes, fill = condition) +
      geom_boxplot(alpha = 0.3) +
      geom_point() +
      scale_fill_brewer(palette = "Accent", direction = 1) +
      ggthemes::theme_few() +
      labs(x ="",
           y = "NES") +
      theme(
        legend.position = "top",
        axis.text.x = element_text(angle = 90)
      ) +
      facet_wrap(vars(process), scales = "free")
  })
  
  
  
  
  #Gene Overlap All DEGs
  df_overlap <- reactive({
    req(input$go_overlap)
    req(conditions_genes())
    
    conditions_genes() %>%
      filter(log2fold_change >= input$filter_logfc_overlap | log2fold_change <= - input$filter_logfc_overlap,
             padj < input$filter_padj_overlap) %>% 
      dplyr::select(process = condition, genes)
  })
  
  # Reactive expression to compute GO enrichment results
  overlap_results <- reactive({
    req(input$go_overlap)
    req(df_overlap())
    
    # Initialize variables
    unique_cond <- unique(df_overlap()$process) # List of unique conditions
    shared_genes_df <- data.frame() # Data frame to store shared gene results
    
    # Compute overlap between conditions
    for (i in 1:(length(unique_cond) - 1)) {
      for (j in (i + 1):length(unique_cond)) {
        resultado_temp <- overlap_genes(unique_cond[i], unique_cond[j], df_overlap())
        shared_genes_df <- bind_rows(shared_genes_df, resultado_temp)
      }
    }
    
    # Compute Fisher's exact test p-values
    shared_genes_df <- shared_genes_df %>%
      mutate(pvalue = pmap_dbl(list(Shared, NotShared_cond1, NotShared_cond2), fisher_exact_test)) %>% 
      mutate("-log(p)" = -log10(pvalue))
    
    # Create mirror table
    shared_genes_df_mirror <- shared_genes_df %>% 
      rename(Cond1_temp = Cond1, Cond2_temp = Cond2) %>%
      rename(Cond1 = Cond2_temp, Cond2 = Cond1_temp) %>%
      bind_rows(shared_genes_df) %>%
      mutate(Percentage_Shared_Cond1 = round(Percentage_Shared_Cond1, 2),
             Percentage_Shared_Cond2 = round(Percentage_Shared_Cond2, 2))
    
    # Calculate Jaccard distance matrix
    jaccard.matrix <- df_overlap() %>%
      mutate(present = 1) %>% 
      distinct() %>% 
      pivot_wider(names_from = genes, values_from = present, values_fill = 0) %>%
      column_to_rownames(var = "process") %>%
      vegdist(binary = TRUE, method = "jaccard", diag = TRUE, upper = TRUE) %>%
      as.matrix() %>%
      {1 - .}
    
    # Combine Jaccard distance with shared genes data
    shared_genes_df_mirror <- jaccard.matrix %>% 
      as.data.frame() %>% 
      rownames_to_column("Cond1") %>% 
      pivot_longer(cols = -"Cond1", names_to = "Cond2", values_to = "jaccard_distance") %>% 
      inner_join(shared_genes_df_mirror, by = c("Cond1", "Cond2")) %>% 
      mutate(jaccard_distance = round(jaccard_distance, 3))
    
    return(shared_genes_df_mirror)
  })
  
  # Render the final results in a DataTable
  output$overlap_table <- DT::renderDataTable({
    req(input$go_overlap)
    req(overlap_results())
    
    
    overlap_results() %>% 
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
  
  #Gene Overlap All DEGs Plot
  output$overlap_heatmap <- renderPlotly({
    req(overlap_results())
    req(input$overlap_metric) # Garante que a escolha da métrica está disponível
    
    # Prepare the data
    heatmap_data <- overlap_results()
    
    # Ensure that the data is in a matrix format required by heatmaply
    heatmap_matrix <- heatmap_data %>%
      dplyr::select(Cond1, Cond2, overlap_metric = input$overlap_metric) %>%
      pivot_wider(names_from = Cond2, values_from = overlap_metric) %>%
      distinct() %>%
      column_to_rownames(var = "Cond1") %>%
      as.matrix()
    
    # Create the heatmap
    heatmaply(heatmap_matrix,
              xlab = "Conditions",
              ylab = "Conditions",
              main = "Heatmap of All Overlapped Genes",
              plot_method = "plotly",
              colors = c("white", "#FF6666")) %>%
      layout(title = "Heatmap of All Overlapped Genes",
             xaxis = list(title = "Condition"),
             yaxis = list(title = "Condition"))
  })
  
  #Gene Overlap Immune
  df_immune <- reactive({
    req(input$go_overlap)
    req(conditions_genes())
    
    conditions_genes() %>%
      filter(log2fold_change >= input$filter_logfc_overlap | 
               log2fold_change <= - input$filter_logfc_overlap,
             padj < input$filter_padj_overlap) %>% 
      dplyr::select(process = condition, genes) %>% 
      inner_join(ImmuneGO_genes_general %>% dplyr::select(genes) %>% distinct(), 
                 by = "genes")
  })
  
  # Reactive expression to compute GO enrichment results
  overlap_results_immune <- reactive({
    req(input$go_overlap)
    req(df_immune())
    
    # Initialize variables
    unique_cond <- unique(df_overlap()$process) # List of unique conditions
    shared_genes_df <- data.frame() # Data frame to store shared gene results
    
    # Compute overlap between conditions
    for (i in 1:(length(unique_cond) - 1)) {
      for (j in (i + 1):length(unique_cond)) {
        resultado_temp <- overlap_genes(unique_cond[i], unique_cond[j], df_immune())
        shared_genes_df <- bind_rows(shared_genes_df, resultado_temp)
      }
    }
    
    # Compute Fisher's exact test p-values
    shared_genes_df <- shared_genes_df %>%
      mutate(pvalue = pmap_dbl(list(Shared, NotShared_cond1, NotShared_cond2), fisher_exact_test)) %>% 
      mutate("-log(p)" = -log10(pvalue))
    
    # Create mirror table
    shared_genes_df_mirror <- shared_genes_df %>% 
      rename(Cond1_temp = Cond1, Cond2_temp = Cond2) %>%
      rename(Cond1 = Cond2_temp, Cond2 = Cond1_temp) %>%
      bind_rows(shared_genes_df) %>%
      mutate(Percentage_Shared_Cond1 = round(Percentage_Shared_Cond1, 2),
             Percentage_Shared_Cond2 = round(Percentage_Shared_Cond2, 2))
    
    # Calculate Jaccard distance matrix
    jaccard.matrix <- df_overlap() %>%
      mutate(present = 1) %>% 
      distinct() %>% 
      pivot_wider(names_from = genes, values_from = present, values_fill = 0) %>%
      column_to_rownames(var = "process") %>%
      vegdist(binary = TRUE, method = "jaccard", diag = TRUE, upper = TRUE) %>%
      as.matrix() %>%
      {1 - .}
    
    # Combine Jaccard distance with shared genes data
    shared_genes_df_mirror <- jaccard.matrix %>% 
      as.data.frame() %>% 
      rownames_to_column("Cond1") %>% 
      pivot_longer(cols = -"Cond1", names_to = "Cond2", values_to = "jaccard_distance") %>% 
      inner_join(shared_genes_df_mirror, by = c("Cond1", "Cond2")) %>% 
      mutate(jaccard_distance = round(jaccard_distance, 3))
    
    return(shared_genes_df_mirror)
  })
  
  # Render the final results in a DataTable
  output$overlap_table_immune <- DT::renderDataTable({
    req(input$go_overlap)
    req(overlap_results_immune())
    
    overlap_results_immune() %>% 
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
  
  #Gene Overlap Immune Plot
  output$overlap_heatmap_immune <- renderPlotly({
    req(overlap_results_immune())
    req(input$overlap_metric) # Garante que a escolha da métrica está disponível
    
    # Prepare the data
    heatmap_data <- overlap_results_immune()
    
    # Ensure that the data is in a matrix format required by heatmaply
    heatmap_matrix <- heatmap_data %>%
      dplyr::select(Cond1, Cond2, overlap_metric = input$overlap_metric) %>%
      pivot_wider(names_from = Cond2, values_from = overlap_metric) %>%
      distinct() %>%
      column_to_rownames(var = "Cond1") %>%
      as.matrix()
    
    # Create the heatmap
    heatmaply(heatmap_matrix,
              xlab = "Conditions",
              ylab = "Conditions",
              main = "Heatmap of Overlapped immune-related Genes",
              plot_method = "plotly",
              colors = c("white", "#FF6666")) %>%
      layout(title = "Heatmap of Overlapped immune-related Genes",
             xaxis = list(title = "Condition"),
             yaxis = list(title = "Condition"))
  })
  
  #Gene Overlap Not Immune
  df_not_immune <- reactive({
    req(input$go_overlap)
    req(conditions_genes())
    
    conditions_genes() %>%
      filter(log2fold_change >= input$filter_logfc_overlap | 
               log2fold_change <= - input$filter_logfc_overlap,
             padj < input$filter_padj_overlap) %>% 
      dplyr::select(process = condition, genes) %>% 
      anti_join(ImmuneGO_genes_general %>% dplyr::select(genes) %>% distinct(), 
                by = "genes")
  })
  
  # Reactive expression to compute GO enrichment results
  overlap_results_not_immune <- reactive({
    req(input$go_overlap)
    req(df_not_immune())
    
    # Initialize variables
    unique_cond <- unique(df_overlap()$process) # List of unique conditions
    shared_genes_df <- data.frame() # Data frame to store shared gene results
    
    # Compute overlap between conditions
    for (i in 1:(length(unique_cond) - 1)) {
      for (j in (i + 1):length(unique_cond)) {
        resultado_temp <- overlap_genes(unique_cond[i], unique_cond[j], df_not_immune())
        shared_genes_df <- bind_rows(shared_genes_df, resultado_temp)
      }
    }
    
    # Compute Fisher's exact test p-values
    shared_genes_df <- shared_genes_df %>%
      mutate(pvalue = pmap_dbl(list(Shared, NotShared_cond1, NotShared_cond2), fisher_exact_test)) %>% 
      mutate("-log(p)" = -log10(pvalue))
    
    # Create mirror table
    shared_genes_df_mirror <- shared_genes_df %>% 
      rename(Cond1_temp = Cond1, Cond2_temp = Cond2) %>%
      rename(Cond1 = Cond2_temp, Cond2 = Cond1_temp) %>%
      bind_rows(shared_genes_df) %>%
      mutate(Percentage_Shared_Cond1 = round(Percentage_Shared_Cond1, 2),
             Percentage_Shared_Cond2 = round(Percentage_Shared_Cond2, 2))
    
    # Calculate Jaccard distance matrix
    jaccard.matrix <- df_overlap() %>%
      mutate(present = 1) %>% 
      distinct() %>% 
      pivot_wider(names_from = genes, values_from = present, values_fill = 0) %>%
      column_to_rownames(var = "process") %>%
      vegdist(binary = TRUE, method = "jaccard", diag = TRUE, upper = TRUE) %>%
      as.matrix() %>%
      {1 - .}
    
    # Combine Jaccard distance with shared genes data
    shared_genes_df_mirror <- jaccard.matrix %>% 
      as.data.frame() %>% 
      rownames_to_column("Cond1") %>% 
      pivot_longer(cols = -"Cond1", names_to = "Cond2", values_to = "jaccard_distance") %>% 
      inner_join(shared_genes_df_mirror, by = c("Cond1", "Cond2")) %>% 
      mutate(jaccard_distance = round(jaccard_distance, 3))
    
    return(shared_genes_df_mirror)
  })
  
  # Render the final results in a DataTable
  output$overlap_table_not_immune <- DT::renderDataTable({
    req(input$go_overlap)
    req(overlap_results_not_immune())
    
    overlap_results_not_immune() %>% 
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
  
  #Gene Overlap Not Immune Plot
  output$overlap_heatmap_not_immune <- renderPlotly({
    req(overlap_results_not_immune())
    req(input$overlap_metric) # Garante que a escolha da métrica está disponível
    
    # Prepare the data
    heatmap_data <- overlap_results_not_immune()
    
    # Ensure that the data is in a matrix format required by heatmaply
    heatmap_matrix <- heatmap_data %>%
      dplyr::select(Cond1, Cond2, overlap_metric = input$overlap_metric) %>%
      pivot_wider(names_from = Cond2, values_from = overlap_metric) %>%
      distinct() %>%
      column_to_rownames(var = "Cond1") %>%
      as.matrix()
    
    # Create the heatmap
    heatmaply(heatmap_matrix,
              xlab = "Conditions",
              ylab = "Conditions",
              main = "Heatmap of Overlapped non-immune Genes",
              plot_method = "plotly",
              colors = c("white", "#FF6666")) %>%
      layout(title = "Heatmap of Overlapped non-immune Genes",
             xaxis = list(title = "Condition"),
             yaxis = list(title = "Condition"))
  })
  
}

shinyApp(ui, server)