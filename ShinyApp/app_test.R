#Import library ----
# Core Shiny and data manipulation libraries
library(shiny)
library(Matrix)
library(tidyverse) # Includes dplyr, ggplot2, tidyr, etc. for data handling
library(janitor)   # For cleaning data names (e.g., clean_names())
library(readxl)    # For reading Excel files
library(scales)    # For custom scale manipulation

# Plotting and visualization
library(plotly)    # For interactive plots (Volcano, Bar, Heatmaps)
library(ggsci)     # Scientific journal color palettes for ggplot2
library(patchwork) # For combining multiple ggplot2 plots
library(ggthemes)  # Extra ggplot2 themes
library(circlize)   # For circular visualization (not explicitly used here, but loaded)
library(RColorBrewer)# For color palettes
library(heatmaply) # For interactive heatmaps

# Bioconductor/Bioinformatics libraries
library(org.Hs.eg.db) # Annotation package for Homo sapiens (used for GSEA/ORA)
library(clusterProfiler) # Core package for GSEA and ORA
library(corto)       # Likely used for ssGSEA (single-sample GSEA)
library(BiocManager) # Bioconductor package management
library(devtools)    # R package development tools (often used for installing from GitHub)

# Statistical and advanced analysis
library(vegan)     # Community ecology package (used here for Jaccard distance calculation in overlap)
library(DT)        # DataTables for interactive R tables
library(BH)        # Boost C++ header files (dependency for some R packages)
library(cytolib)   # Flow cytometry data structures (dependency, but not explicitly used here)

# Project organization
library(here)      # For robust file path specification

#Increase upload file size
options(shiny.maxRequestSize=200*1024^2) # Sets max file upload size to 200MB

################### Custom Functions -------
##### autoORA (for multiple conditions)
# Function to perform Over-Representation Analysis (ORA) across multiple experimental conditions.
autoORA <- function(df,
                    TERM2GENE, # A data frame of gene set terms and their corresponding genes
                    pAdjustMethod, 
                    pvalueCutoff_ora, 
                    qvalueCutoff_ora,
                    geneset_name) {
  
  # Inner function to perform enrichment for a single condition
  enricher_conditional <- function(condition_2, df, TERM2GENE, pAdjustMethod, pvalueCutoff_ora, qvalueCutoff_ora) {
    # Extract genes for the current condition
    genes <- df %>% 
      filter(condition == condition_2) %>% 
      pull(genes)
    
    # Use tryCatch to handle errors (e.g., if a condition has no genes) gracefully
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
      NULL # Returns NULL if an error occurs
    })
  }
  
  # Apply the enrichment function to all unique conditions in the input data frame
  ora_results <- df %>%
    distinct(condition) %>%
    pull(condition) %>%
    map_dfr(~ enricher_conditional(.x, df, TERM2GENE, pAdjustMethod, pvalueCutoff_ora, qvalueCutoff_ora))
  
  # Final formatting of the results
  ora_results_final <- ora_results %>%
    rownames_to_column("delete") %>% # Removes rownames which are IDs from enricher
    dplyr::select(-delete, -Description) %>%
    dplyr::select(condition, geneset_name, everything()) %>% 
    clean_names() %>% # Converts column names to snake_case (e.g., pvalueCutoff -> pvalue_cutoff)
    rename(process = id) # Renames the ID column to process
  
  return(ora_results_final)
}


####### autoGSEA (for multiple conditions)
# Function to perform Gene Set Enrichment Analysis (GSEA) across multiple experimental conditions.
autoGSEA <- function(df, TERM2GENE, 
                     geneset_name, 
                     minGSSize, 
                     maxGSSize,
                     pvalueCutoff,
                     pAdjustMethod) {
  results <- list()
  
  # Get all unique conditions
  condicoes <- df %>%
    dplyr::pull(condition) %>%
    unique() %>%
    as.character()
  
  # Loop through each condition
  for (condicao in condicoes) {
    # Prepare the ranked gene list for GSEA for the current condition
    degs_condicao <- df %>%
      dplyr::filter(condition == condicao) %>%
      dplyr::select(genes, log2fold_change) %>%
      distinct() %>%
      # Rank by log2fold_change (ties are broken randomly)
      mutate(rank = rank(log2fold_change, ties.method = "random")) %>%
      arrange(desc(rank))
    
    # Create the named numeric vector required by GSEA (log2FC as values, gene names as names)
    gene_list_lf2c <- as.vector(degs_condicao$log2fold_change)
    names(gene_list_lf2c) <- degs_condicao$genes
    gene_list_lf2c <- na.omit(gene_list_lf2c)
    gene_list_lf2c <- sort(gene_list_lf2c, decreasing = TRUE)
    
    # Perform GSEA, handling potential errors
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
    }, 
    error = function(e) {
      print(paste("Error in GSEA for condition:", condicao, "Error:", e$message))
      return(NULL) # Returns NULL on error
    })
    
    # Store results if successful
    if (!is.null(auto_gsea)) {
      results[[paste(condicao, geneset_name, sep = "_")]] <- auto_gsea
    }
  }
  
  # Combine and format all results
  final_result <- bind_rows(results) %>%
    rename(process = ID) %>% # Rename ID to process
    dplyr::select(-Description) %>%
    clean_names() %>%
    dplyr::select(condition, everything())
  
  return(final_result)
}


# Function for overlapping -------
# Calculates the overlap between two gene sets and common overlap metrics.
overlap_genes <- function(cond1, cond2, data) {
  # Extract genes for condition 1 and condition 2
  genes_cond1 <- data$genes[data$process == cond1]
  genes_cond2 <- data$genes[data$process == cond2]
  
  # Calculate intersection (shared genes)
  genes_shared <- intersect(genes_cond1, genes_cond2)
  
  # Calculate unique genes for each condition
  genes_notshared_cond1 <- setdiff(genes_cond1, genes_cond2)
  genes_notshared_cond2 <- setdiff(genes_cond2, genes_cond1)
  
  # Total gene counts
  total_genes_cond1 <- length(genes_cond1)
  total_genes_cond2 <- length(genes_cond2)
  
  # Calculate percentage shared
  percentage_shared_cond1 <- length(genes_shared) / total_genes_cond1 * 100
  percentage_shared_cond2 <- length(genes_shared) / total_genes_cond2 * 100
  
  # Compile results into a data frame
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
# Performs Fisher's exact test on the 2x2 contingency table for gene overlap.
fisher_exact_test <- function(shared, notshared1, notshared2) {
  # Contingency table:
  # | Shared | Not Shared 1 |
  # | Not Shared 2 | Other/Universe |
  # The 'Other/Universe' cell is set to 0 here as the background universe size is not passed to this function,
  # simplifying it for relative comparison between overlaps (a common pattern in overlap analysis).
  cont_table <- matrix(c(shared, notshared1, notshared2, 0), nrow = 2)
  results_fisher <- fisher.test(cont_table)
  return(results_fisher$p.value)
}

# Load data ------
# This section pre-loads immunological gene sets and annotation data from GitHub.

# Gene sets 

######## IMMUNEGO 
# Load ImmuneGO Annotated Gene Ontology data
ImmuneGO_Annotated_GO <- readRDS(url("https://github.com/wapsyed/VaxGO/raw/main/Tables/ImmuneGO_Annotated_GO_2024-05-28.rds", method="libcurl"))

# Annotation table for ImmuneGO
ImmuneGO_annotation = ImmuneGO_Annotated_GO %>% 
  dplyr::select(process, immune_system:immune_tissue) 

# General ImmuneGO gene sets (manual terms, excluding broad terms)
ImmuneGO_genes_general = ImmuneGO_Annotated_GO %>% 
  separate_rows(genes, sep = ",") %>% 
  filter(go_term == "Manual",
         !process %in% c("ADAPTIVE IMMUNE SYSTEM",
                         "INNATE IMMUNE SYSTEM",
                         "BCR REPERTOIRE",
                         "TCR REPERTOIRE")) %>% 
  dplyr::select(process, genes)

# Specific ImmuneGO gene sets (non-manual GO terms)
ImmuneGO_genes_specific = ImmuneGO_Annotated_GO %>% 
  separate_rows(genes, sep = ",") %>% 
  filter(!go_term == "Manual") %>% 
  dplyr::select(process=gene_set_short, genes)

######## CELL MARKER IMMUNE 
# Load CellMarker Immune Cell gene sets
CellMarker_ImmuneCells = readRDS(url("https://github.com/wapsyed/VaxGO/raw/main/Tables/CellMarker_ImmuneCells.rds", method="libcurl")) %>% 
  dplyr::rename(process = cell_name, genes = marker)

# Annotation table for CellMarker
CellMarker_annotation = CellMarker_ImmuneCells %>% 
  dplyr::select(process, Type) %>% 
  distinct()

# Gene sets for CellMarker
CellMarker_genes = CellMarker_ImmuneCells %>% 
  dplyr::select(-Type)

######## VAX MSIGDB 
# Load Vaccine Signature Database (VaxSigDB) gene sets
VaxSigDB_Gene_sets_Annotated_RAW <- readRDS(url("https://github.com/wapsyed/VaxGO/raw/main/Tables/VaxSigDB_Gene_sets_Annotated_RAW.rds", method = "libcurl"))

# Filters to be applied
vaxsig_columns = c("VACCINE", "PLATFORM", "TARGET_PATHOGEN_DISEASE", "MICROBE_TYPE", "STUDY_TYPE", "STUDY_SUBTYPE","DATE-TIME", "DATE", "TIME", "AGE", "AGE_CATEGORY", "SAMPLE_SOURCE", "SYSTEMATIC_NAME", "GENE_SYMBOLS")
filter_vaxsig_samplesource = "PBMC" # Filter for Peripheral Blood Mononuclear Cells
filter_vaxsig_studytype = "VACCINE"
filter_vaxsig_studysubtype = "VAC ONLY"

# Filtered VaxSigDB gene sets
VaxSigDB_Genesets_filtered = VaxSigDB_Gene_sets_Annotated_RAW %>% 
  dplyr::select(vaxsig_columns) %>% 
  filter(SAMPLE_SOURCE == filter_vaxsig_samplesource,
         STUDY_TYPE == filter_vaxsig_studytype,
         STUDY_SUBTYPE == filter_vaxsig_studysubtype) 
# Annotation for filtered VaxSigDB
VaxSigDB_annotation = VaxSigDB_Genesets_filtered %>% 
  mutate(process = paste0(VACCINE, " (", `DATE-TIME`, ", ", AGE, " YO)")) %>% 
  dplyr::select(process, PLATFORM, TARGET_PATHOGEN_DISEASE, MICROBE_TYPE, DATE, TIME)
# Genes for filtered VaxSigDB
VaxSigDB_Genes_filtered = VaxSigDB_Genesets_filtered %>% 
  separate_rows(GENE_SYMBOLS, sep = ",") %>% 
  mutate(process = paste0(VACCINE, " (", `DATE-TIME`, ", ", AGE, " YO)")) %>% 
  dplyr::select(process, genes = GENE_SYMBOLS)

######## BTM MODULES 
# Load Blood Transcriptional Modules (BTM) annotation table
btm_annotation_table <- readRDS(url("https://github.com/wapsyed/VaxGO/raw/main/Tables/btm_annotation_immune_table_labelled.rds", method = "libcurl"))

# Filters to be applied
filter_btm_category = "immune"
filter_btm_annotationlevel = c("complete", "partial") # Filter for 'immune' category and 'complete' or 'partial' annotation
# filter_btm_annotationlevel = "complete" # Commented out alternative filter
# filter_btm_annotationlevel = "partial" # Commented out alternative filter

# Annotation for filtered BTM
btm_annotation = btm_annotation_table %>% 
  filter(module_category == filter_btm_category,
         annotation_level %in% filter_btm_annotationlevel) %>% 
  dplyr::select(immune_system, immune_subsystem, cells, composite_name)
# Genes for filtered BTM
btm_genes = btm_annotation_table %>% 
  filter(module_category == filter_btm_category,
         annotation_level %in% filter_btm_annotationlevel) %>% 
  dplyr::select(process = composite_name, genes = module_member_genes) %>% 
  separate_rows(genes, sep = ",")

# Vaccine Atlas (Example DEGs)
# Load a pre-filtered example dataset of Differentially Expressed Genes (DEGs)
# Only DEGs, p.adj <0.05, no Log2FC filter.
covid_13vax = readRDS(url("https://github.com/wapsyed/VaxGO/raw/main/Tables/degs_covid_vaxhagan.rds", method = "libcurl"))

# Filtered example DEG data
covid_13vax_Genes_filtered = covid_13vax %>% 
  select(condition, genes)



# UI definition ------
ui <- fluidPage(
  # Custom CSS styles for aesthetics and font loading
  tags$head(
    tags$style(HTML("
      body {
        background-color: #ffffff;
        color: #000000;
        font-family: 'Montserrat', sans-serif; /* Set primary font */
        font-size: 12px;
      }
      .form-group,
      .shiny-input-container {
        font-size: 12px;
      }
      /* Ensure input elements also use the correct font size */
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
      /* Adjust font size for tabs */
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
    # Link to load the Montserrat font
    tags$link(
      href = "https://fonts.googleapis.com/css2?family=Montserrat:wght@400;700&display=swap",
      rel = "stylesheet"
    ),
    # Links to load Font Awesome and Ionicons for social media icons
    tags$link(
      rel = "stylesheet",
      href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css"
    ),
    tags$link(
      rel = "stylesheet",
      href = "https://code.ionicframework.com/ionicons/2.0.1/css/ionicons.min.css"
    )
  ),
  # Header containing the title and social media links
  div(
    style = "display: flex; justify-content: space-between; align-items: center;",
    titlePanel(
      div(
        h2("VaxGO Tool", style = "margin-bottom: 0px;"),
        # Developer attribution and GitHub link
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
    # Social media navigation links with icons
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
  # Main layout: Sidebar and Main Panel
  sidebarLayout(
    # Sidebar Panel (3 units wide) for all inputs/parameters
    sidebarPanel(
      width = 3,
      h4("Parameters"), # Title for the sidebar
      # Tabset for organizing input parameters
      tabsetPanel(
        type = "tabs",
        # DEGs Analysis Parameters Tab
        tabPanel("DEGs",
                 # Link to download example file
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
                 # Input for uploading the main DEG file
                 fileInput("conditions_genes", 
                           "Upload your file",
                           accept = ".csv"),
                 
                 
                 h4("Filters", style = "margin-top: 20px;"),
                 # Input for adjusted p-value cutoff
                 numericInput(inputId = "filter_padj", 
                              label = "Adjusted p-value cutoff", 
                              value = 0.05, 
                              min = 0, max = 1, step = 0.01),
                 
                 # Input for Log2 Fold Change cutoff
                 numericInput(inputId = "filter_logfc", 
                              label = "Log2 Fold Change cutoff", 
                              value = 1, 
                              min = 0, max = 10, step = 0.1),
                 # Action button to trigger initial DEG processing/plotting
                 actionButton("go_1", "Go")
        ),
        # ORA (Over-Representation Analysis) Parameters Tab
        tabPanel("ORA",
                 h4("ORA Parameters", style = "margin-top: 20px;"),
                 
                 # Input for ORA p-value cutoff
                 numericInput(inputId = "pvalueCutoff_ora", 
                              label = "ORA p-value Cutoff", 
                              value = 0.25, 
                              min = 0, max = 1, step = 0.01),
                 
                 # Input for ORA q-value cutoff
                 numericInput(inputId = "qvalueCutoff_ora", 
                              label = "ORA q-value Cutoff", 
                              value = 0.25, 
                              min = 0, max = 1, step = 0.01),
                 
                 # Dropdown for p-value adjustment method
                 selectInput(inputId = "pAdjustMethod_ora", 
                             label = "p-value Adjustment Method", 
                             choices = c("BH", "BY", "fdr", "holm", "hochberg", "hommel", "bonferroni"),
                             selected = "BH"),
                 
                 # Action button to trigger ORA
                 actionButton("go_ora", "Go")
                 
        ),
        # GSEA (Gene Set Enrichment Analysis) Parameters Tab
        tabPanel("GSEA",
                 h4("GSEA Parameters", style = "margin-top: 20px;"),
                 
                 # Input for minimum gene set size
                 numericInput(inputId = "minGSSize", 
                              label = "Minimum Gene Set Size", 
                              value = 1, 
                              min = 1, max = 1000, step = 1),
                 
                 # Input for maximum gene set size
                 numericInput(inputId = "maxGSSize", 
                              label = "Maximum Gene Set Size", 
                              value = 1000, 
                              min = 1, max = 10000, step = 1),
                 
                 # Dropdown for p-value adjustment method
                 selectInput(inputId = "pAdjustMethod", 
                             label = "p-value Adjustment Method", 
                             choices = c("BH", "BY", "fdr", "holm", "hochberg", "hommel", "bonferroni"),
                             selected = "BH"),
                 
                 # Input for GSEA p-value cutoff
                 numericInput(inputId = "pvalueCutoff_gsea", 
                              label = "GSEA p-value Cutoff", 
                              value = 0.25, 
                              min = 0, max = 1, step = 0.01),
                 
                 # Dropdown for organism database
                 selectInput(inputId = "organism", 
                             label = "Organism", 
                             choices = c("org.Hs.eg.db", "org.Mm.eg.db"),
                             selected = "org.Hs.eg.db"),
                 
                 # Action button to trigger GSEA
                 actionButton("go_2", "Go")
        ),
        # ssGSEA (Single-Sample GSEA) Parameters Tab
        tabPanel("ssGSEA",
                 h4("ssGSEA Parameters", style = "margin-top: 20px;"),
                 
                 # Link to download example normalized counts file
                 tags$p(
                   HTML(
                     paste0(
                       "Download an ",
                       "<a href='https://github.com/wapsyed/VaxGO/blob/662205136a565e3420c26bfce665073641c45cc4/Example/Example_normalized_counts.csv'>",
                       "example"
                     )
                   ),
                 ),
                 # Input for normalized counts data
                 fileInput(inputId = "ssgsea_gene", 
                           label = "Normalized counts",
                           accept = ".csv"),
                 
                 
                 # Link to download example sample annotation file (metadata)
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
                 # Input for sample metadata
                 fileInput(inputId = "metadata", 
                           label = "Sample annotations",
                           accept = ".csv"),
                 
                 # Action button to trigger ssGSEA calculation
                 actionButton("go_ssgsea", "Go"),
                 # Dropdown to select a cell/process for plotting ssGSEA results
                 selectInput("selected_cellmarker", 
                             "Select a cell:", 
                             choices = NULL,  # Initially empty, updated dynamically
                             selected = NULL,
                             multiple = FALSE)
        ),
        # Gene Overlap Parameters Tab
        tabPanel("Overlap",
                 h4("Overlap Parameters", style = "margin-top: 20px;"),
                 # DEG filter for overlap calculation (adjusted p-value)
                 numericInput(inputId = "filter_padj_overlap", 
                              label = "Adjusted p-value cutoff", 
                              value = 0.05, 
                              min = 0, max = 1, step = 0.01),
                 
                 # DEG filter for overlap calculation (Log2 Fold Change)
                 numericInput(inputId = "filter_logfc_overlap", 
                              label = "Log2 Fold Change cutoff", 
                              value = 1, 
                              min = 0, max = 10, step = 0.1),
                 
                 # Radio buttons to choose the metric for heatmap visualization
                 radioButtons(inputId = "overlap_metric", 
                              label = "Overlap method", 
                              choices = c("jaccard_distance", "Shared"),
                              selected = "jaccard_distance"),
                 
                 # Action button to trigger overlap analysis
                 actionButton("go_overlap", "Go")
        )
      )
    ),
    
    #RESULTS ------
    # Main Panel (9 units wide) for all outputs
    mainPanel(
      width = 9,
      h4("Results"), # Title for the results section
      
      # Tabset for organizing output results
      tabsetPanel(
        type = "tabs",
        # Tab for DEGs Analysis results
        tabPanel("DEGs Analysis",
                 fluidRow(
                   column(width = 6,
                          # Interactive Volcano Plot
                          plotlyOutput("volcano_plot")
                   ),
                   column(width = 6,
                          # Interactive Bar Plot of DEG counts
                          plotlyOutput("bar_plot")
                   )
                 ),
                 br(),
                 h4("Genes Table"),
                 # Interactive table of all significant DEGs
                 DT::dataTableOutput("degs_table")
        ),
        # Tab for ORA results
        tabPanel("ORA",
                 tabsetPanel(
                   type = "tabs",
                   # Sub-tab for ImmuneGO General ORA results
                   tabPanel("ImmuneGO General",
                            DT::dataTableOutput("ora_table_immuneGO_general"),
                            plotlyOutput("ora_plot_immuneGO_general")),
                   tabPanel("ImmuneGO Specific",
                            DT::dataTableOutput("ora_table_immuneGO_specific"),
                            plotlyOutput("ora_plot_immuneGO_specific")),
                   tabPanel("CellMarker Immune",
                            DT::dataTableOutput("gsea_table_cellmarker"),
                            plotlyOutput("ora_plot_cellmarker")),
                   tabPanel("BTM immune",
                            DT::dataTableOutput("gsea_table_btm_immune"),
                            plotlyOutput("ora_plot_btm_immune")),
                   tabPanel("Vax MSigDB",
                            DT::dataTableOutput("gsea_table_vaxsigdb"),
                            plotlyOutput("ora_plot_vaxsigdb"))
                 )
        ),
        # Tab for GSEA results
        tabPanel("GSEA",
                 tabsetPanel(
                   type = "tabs",
                   # Sub-tabs for GSEA results for different gene sets
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
        # Tab for ssGSEA results
        tabPanel("ssGSEA",
                 DT::dataTableOutput("ssgsea_table_cellmarker"),
                 plotlyOutput("ssgsea_plot_cellmarker")
        ),
        # Tab for Gene Overlap results
        tabPanel("Gene Overlap",
                 # Heatmap and table for ALL DEGs overlap
                 plotlyOutput("overlap_heatmap"),
                 DT::dataTableOutput("overlap_table"),
                 br(),
                 # Heatmap and table for IMMUNE-RELATED DEGs overlap
                 plotlyOutput("overlap_heatmap_immune"),
                 DT::dataTableOutput("overlap_table_immune"),
                 br(),
                 # Heatmap and table for NON-IMMUNE DEGs overlap
                 plotlyOutput("overlap_heatmap_not_immune"),
                 DT::dataTableOutput("overlap_table_not_immune"),
        )
      )
    )
  )
)

















# Server definition ------
server <- function(input, output, session) {
  
  # Reactive to load and process the uploaded DEG file ----
  conditions_genes <- reactive({
    req(input$conditions_genes) # Ensure a file has been uploaded
    
    # Load the file
    conditions_genes <- read.csv(input$conditions_genes$datapath) %>% 
      distinct() # Remove duplicate rows
    conditions_genes
  })
  
  # Reactive to filter for significant DEGs based on p-value ----
  degs_sig <- reactive({
    req(conditions_genes(), input$filter_padj) # Ensure data and filter are available
    
    degs_sig <- conditions_genes() %>%
      filter(padj < input$filter_padj) %>% # Filter by adjusted p-value cutoff
      distinct()
    
    degs_sig
  })
  
  # Reactive to filter DEGs by both p-value and Log2FC (for gene lists) -----
  genes <- reactive({
    req(degs_sig())
    
    genes <- degs_sig() %>%
      filter(log2fold_change >= input$filter_logfc | log2fold_change <= - input$filter_logfc, # Filter by Log2FC cutoff
             padj < input$filter_padj) %>% # Filter by adjusted p-value cutoff
      pull(genes) # Extract only the gene names
    
    genes
  })
  
  # Volcano plot output (DEG Analysis Tab) -----
  output$volcano_plot <- renderPlotly({
    req(input$go_1) # Requires the "Go" button in the DEGs tab to be pressed
    req(conditions_genes())
    
    # Prepare the data: calculate -log10(padj)
    data <- conditions_genes() %>%
      mutate(log_q = -log10(padj))
    
    # Get unique conditions for dropdown menu (to switch between conditions)
    unique_conditions <- unique(data$condition)
    
    # Create an empty plot object
    p <- plot_ly()
    
    # Add a trace (data series) for each condition
    for (cond in unique_conditions) {
      p <- p %>%
        add_trace(
          data = data[data$condition == cond, ],
          x = ~log2fold_change, 
          y = ~log_q, 
          text = ~genes, # Gene name on hover
          type = 'scatter', 
          mode = 'markers',
          marker = list(
            color = ~log2fold_change, # Color by Log2FC
            size = 10,
            opacity = 0.9,
            colorscale = c("#3399FF", "white", "#FF6666") # Custom colorscale (blue, white, red)
          ),
          name = cond,
          visible = ifelse(cond == unique_conditions[1], TRUE, FALSE) # Show only the first condition initially
        )
    }
    
    # Create the buttons for switching between conditions
    condition_buttons <- lapply(1:length(unique_conditions), function(i) {
      visibility <- rep(FALSE, length(unique_conditions))
      visibility[i] <- TRUE
      list(
        method = "restyle",
        args = list("visible", visibility),
        label = unique_conditions[i]
      )
    })
    
    # Add layout and the dropdown menu (updatemenus)
    p %>%
      layout(
        title = "Volcano Plot",
        xaxis = list(title = "Log2 Fold Change"),
        yaxis = list(title = "-Log10 Adjusted P-value"),
        updatemenus = list(
          list(
            y = 1,  # Position the buttons near the top
            x = 0,   # Position the buttons horizontally
            xanchor = "center",
            yanchor = "top",
            buttons = condition_buttons,
            showactive = TRUE
          )
        )
      )
  })
  
  # Bar Plot output (DEG Analysis Tab) ------
  output$bar_plot <- renderPlotly({
    req(input$go_1) # Requires the "Go" button in the DEGs tab to be pressed
    req(conditions_genes())
    
    # Prepare the data: filter and count genes by condition and direction (up/downregulated)
    data <- conditions_genes()  %>%
      filter(log2fold_change >= input$filter_logfc | log2fold_change <= - input$filter_logfc,
             padj < input$filter_padj) %>% 
      group_by(condition, direction) %>% 
      summarize(n_direct = n()) %>% # Count DEGs in each direction
      ungroup()
    
    
    # Create the interactive bar plot
    plot_ly(data, 
            x = ~condition, 
            y = ~n_direct, 
            color = ~direction, # Color bars by direction
            colors = c("#3399FF", "#FF6666"), # Blue for down, Red for up
            text = ~direction,
            type = 'bar') %>%
      layout(
        title = "DEGs Counts",
        xaxis = list(title = "Condition"),
        yaxis = list(title = "# Genes"),
        barmode = 'group',  # Group bars side-by-side
        coloraxis = list(colorbar = list(title = 'Direction'))
      )
  })
  
  # DEGs table output (DEG Analysis Tab) -----
  output$degs_table <- DT::renderDataTable({
    req(input$go_1) # Requires the "Go" button in the DEGs tab to be pressed
    req(degs_sig())
    
    # Define input values for the cutoffs for visual highlighting
    filter_l2fc <- input$filter_logfc
    filter_neg_l2fc <- -input$filter_logfc
    
    # Render the interactive DataTable
    degs_sig() %>%
      datatable(options = list(pageLength = 25)) %>%
      # Apply conditional formatting based on Log2 Fold Change
      formatStyle(
        'log2fold_change',
        backgroundColor = styleInterval(
          c(filter_neg_l2fc, filter_l2fc), # Intervals for color transition
          c('#3399FF', 'white', '#FF6666') # Blue (down), white (non-sig), Red (up)
        ),
        fontWeight = 'bold'
      )
  })
  
  
  # ORA ImmuneGO general results ----
  ora_results_ImmuneGO_General <- reactive({
    req(input$go_ora) # Requires the "Go" button in the ORA tab to be pressed
    req(degs_sig())
    
    degs_sig_data <- degs_sig()
    print("Data for ORA ImmuneGO General:") # Logging
    print(head(degs_sig_data))
    
    # Perform ORA using the custom autoORA function
    degs_ORA_ImmuneGO_General <- degs_sig_data %>%
      dplyr::select(condition, genes, log2fold_change) %>% 
      autoORA(TERM2GENE = ImmuneGO_genes_general, # Use general ImmuneGO gene sets
              geneset_name = "ImmuneGO General", 
              pAdjustMethod = input$pAdjustMethod_ora, 
              pvalueCutoff_ora = input$pvalueCutoff_ora,
              qvalueCutoff_ora = input$qvalueCutoff_ora) 
    
    
    print("ORA ImmuneGO General results:") # Logging
    print(head(degs_ORA_ImmuneGO_General))
    
    degs_ORA_ImmuneGO_General
  })
  
  # ORA ImmuneGO general table output 
  output$ora_table_immuneGO_general <- DT::renderDataTable({
    req(input$go_ora)
    req(ora_results_ImmuneGO_General())
    # Render interactive DataTable with Export buttons
    ora_results_ImmuneGO_General() %>% 
      datatable(extensions = 'Buttons',
                options = list(
                  paging = TRUE,
                  searching = TRUE,
                  fixedColumns = TRUE,
                  autoWidth = TRUE,
                  ordering = TRUE,
                  dom = 'Blfrtip', # Layout for buttons, filter, table, pagination, info
                  buttons = c('copy', 'csv', 'excel')
                ),
                
                class = "display"
      )
  })
  
  # ORA ImmuneGO general plot output 
  output$ora_plot_immuneGO_general <- renderPlotly({
    req(input$go_ora) 
    req(ora_results_ImmuneGO_General())
    
    # Prepare the data: filter top 5 enriched processes per condition and reorder
    data <- ora_results_ImmuneGO_General() %>% 
      group_by(condition) %>%
      arrange(condition, desc(-log(qvalue))) %>% 
      slice_head(n = 5) %>% # Select top 5
      ungroup() %>%
      mutate(process = fct_reorder(process, -log(qvalue))) # Reorder processes for plotting
    
    num_conditions <- data %>%
      distinct(condition) %>%
      nrow()
    
    # Create the ggplot bar plot
    plot <- ggplot(data) +
      aes(x = -log(qvalue), y = process,
          fill = -log(qvalue)) +
      geom_col() +
      # Facet by condition, allowing free y-axis for each condition's top 5
      facet_wrap(~condition, scales = "free_y", nrow = num_conditions) +
      labs(x = "-log10(q-value)", y = "Process", title = "ORA Results") +
      theme_minimal()
    
    # Convert ggplot to interactive plotly object
    ggplotly(plot)
  })
  
  # ORA ImmuneGO specific results ----
  ora_results_ImmuneGO_Specific <- reactive({
    req(input$go_ora) # Requires the "Go" button in the ORA tab to be pressed
    req(degs_sig())
    
    degs_sig_data <- degs_sig()
    print("Data for ORA ImmuneGO Specific:") # Logging
    print(head(degs_sig_data))
    
    # Perform ORA using the custom autoORA function
    degs_ORA_ImmuneGO_Specific <- degs_sig_data %>%
      dplyr::select(condition, genes, log2fold_change) %>% 
      autoORA(TERM2GENE = ImmuneGO_genes_specific, # Use general ImmuneGO gene sets
              geneset_name = "ImmuneGO Specific", 
              pAdjustMethod = input$pAdjustMethod_ora, 
              pvalueCutoff_ora = input$pvalueCutoff_ora,
              qvalueCutoff_ora = input$qvalueCutoff_ora) 
    
    
    print("ORA ImmuneGO Specific results:") # Logging
    print(head(degs_ORA_ImmuneGO_Specific))
    
    degs_ORA_ImmuneGO_Specific
  })
  
  # ORA ImmuneGO specific table output 
  output$ora_table_immuneGO_specific <- DT::renderDataTable({
    req(input$go_ora)
    req(ora_results_ImmuneGO_Specific())
    # Render interactive DataTable with Export buttons
    ora_results_ImmuneGO_Specific() %>% 
      datatable(extensions = 'Buttons',
                options = list(
                  paging = TRUE,
                  searching = TRUE,
                  fixedColumns = TRUE,
                  autoWidth = TRUE,
                  ordering = TRUE,
                  dom = 'Blfrtip', # Layout for buttons, filter, table, pagination, info
                  buttons = c('copy', 'csv', 'excel')
                ),
                
                class = "display"
      )
  })
  
  # ORA ImmuneGO specific plot output 
  output$ora_plot_immuneGO_specific <- renderPlotly({
    req(input$go_ora) 
    req(ora_results_ImmuneGO_Specific())
    
    # Prepare the data: filter top 5 enriched processes per condition and reorder
    data <- ora_results_ImmuneGO_Specific() %>% 
      group_by(condition) %>%
      arrange(condition, desc(-log(qvalue))) %>% 
      slice_head(n = 5) %>% # Select top 5
      ungroup() %>%
      mutate(process = fct_reorder(process, -log(qvalue))) # Reorder processes for plotting
    
    num_conditions <- data %>%
      distinct(condition) %>%
      nrow()
    
    # Create the ggplot bar plot
    plot <- ggplot(data) +
      aes(x = -log(qvalue), y = process,
          fill = -log(qvalue)) +
      geom_col() +
      # Facet by condition, allowing free y-axis for each condition's top 5
      facet_wrap(~condition, scales = "free_y", nrow = num_conditions) +
      labs(x = "-log10(q-value)", y = "Process", title = "ORA Results") +
      theme_minimal()
    
    # Convert ggplot to interactive plotly object
    ggplotly(plot)
  })
  
  
  
  
  
  
  
  
  
  
  
  
  
  # GSEA ImmuneGO general results reactive -----
  gsea_results_ImmuneGO_General <- reactive({
    req(input$go_2) # Requires the "Go" button in the GSEA tab to be pressed
    req(degs_sig())
    
    degs_sig_data <- degs_sig()
    print("Data for GSEA ImmuneGO General:") # Logging
    print(head(degs_sig_data))
    
    # Perform GSEA using the custom autoGSEA function
    degs_GSEA_ImmuneGO_General <- degs_sig_data %>%
      dplyr::select(condition, genes, log2fold_change) %>% 
      autoGSEA(TERM2GENE = ImmuneGO_genes_general, # Use general ImmuneGO gene sets
               geneset_name = "ImmuneGO General",
               minGSSize = input$minGSSize,
               maxGSSize = input$maxGSSize,
               pvalueCutoff = input$pvalueCutoff_gsea,
               pAdjustMethod = input$pAdjustMethod) 
    
    print("GSEA ImmuneGO General results:") # Logging
    print(head(degs_GSEA_ImmuneGO_General))
    
    degs_GSEA_ImmuneGO_General
  })
  
  # GSEA ImmuneGO general table output ----
  output$gsea_table_immuneGO_general <- DT::renderDataTable({
    req(input$go_2)
    req(gsea_results_ImmuneGO_General())
    # Render interactive DataTable with Export buttons
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
  
  # GSEA ImmuneGO general plot output ----
  output$gsea_plot_immuneGO_general <- renderPlotly({
    req(input$go_2) 
    req(gsea_results_ImmuneGO_General())
    
    # Prepare the data: filter top 5 enriched processes per condition and reorder
    data <- gsea_results_ImmuneGO_General() %>% 
      group_by(condition) %>%
      arrange(condition, desc(-log(qvalue))) %>% 
      slice_head(n = 5) %>%
      ungroup() %>%
      mutate(process = fct_reorder(process, -log(qvalue)))  
    
    num_conditions <- data %>%
      distinct(condition) %>%
      nrow()
    
    # Create the ggplot bar plot
    plot <- ggplot(data) +
      aes(x = -log(qvalue), y = process,
          fill = -log(qvalue)) +
      geom_col() +
      facet_wrap(~condition, scales = "free_y", nrow = num_conditions) +
      labs(x = "-log10(q-value)", y = "Process", title = "GSEA Results") +
      theme_minimal()
    
    # Convert ggplot to interactive plotly object
    ggplotly(plot)
  })
  
  
  
  # GSEA ImmuneGO specific results reactive ------
  gsea_results_ImmuneGO_genes_specific <- reactive({
    req(input$go_2)
    req(degs_sig())
    
    degs_sig_data <- degs_sig()
    print("Data for GSEA ImmuneGO Specific:")
    print(head(degs_sig_data))
    
    # Perform GSEA using the custom autoGSEA function
    degs_GSEA_ImmuneGO_genes_specific <- degs_sig_data %>%
      dplyr::select(condition, genes, log2fold_change) %>% 
      autoGSEA(TERM2GENE = ImmuneGO_genes_specific, # Use specific ImmuneGO gene sets
               geneset_name = "ImmuneGO Specific",
               minGSSize = input$minGSSize,
               maxGSSize = input$maxGSSize,
               pvalueCutoff = input$pvalueCutoff_gsea,
               pAdjustMethod = input$pAdjustMethod) 
    
    print("GSEA ImmuneGO Specific results:")
    print(head(degs_GSEA_ImmuneGO_genes_specific))
    
    degs_GSEA_ImmuneGO_genes_specific
  })
  
  # GSEA ImmuneGO specific table output ----
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
  
  # GSEA ImmuneGO specific plot output ----
  output$gsea_plot_immuneGO_specific <- renderPlotly({
    req(input$go_2) 
    req(gsea_results_ImmuneGO_genes_specific())
    
    # Prepare the data: filter top 5 enriched processes per condition and reorder
    data <- gsea_results_ImmuneGO_genes_specific()%>% 
      group_by(condition) %>%
      arrange(condition, desc(-log(qvalue))) %>% 
      slice_head(n = 5) %>%
      ungroup() %>%
      mutate(process = fct_reorder(process, -log(qvalue)))  
    
    num_conditions <- data %>%
      distinct(condition) %>%
      nrow()
    
    # Create the ggplot bar plot
    plot <- ggplot(data) +
      aes(x = -log(qvalue), y = process,
          fill = -log(qvalue)) +
      geom_col() +
      facet_wrap(~condition, scales = "free_y", nrow = num_conditions) +
      labs(x = "-log10(q-value)", y = "Process", title = "GSEA Results") +
      theme_minimal()
    
    # Convert ggplot to interactive plotly object
    ggplotly(plot)
  })
  
  # GSEA CellMarker results reactive -----
  gsea_results_CellMarker <- reactive({
    req(input$go_2)
    req(degs_sig())
    
    degs_sig_data <- degs_sig()
    print("Data for GSEA CellMarker:")
    print(head(degs_sig_data))
    
    # Perform GSEA using the custom autoGSEA function
    degs_GSEA_CellMarker <- degs_sig_data %>%
      dplyr::select(condition, genes, log2fold_change) %>% 
      autoGSEA(TERM2GENE = CellMarker_genes, # Use CellMarker gene sets
               geneset_name = "CellMarker immune",
               minGSSize = input$minGSSize,
               maxGSSize = input$maxGSSize,
               pvalueCutoff = input$pvalueCutoff_gsea,
               pAdjustMethod = input$pAdjustMethod) 
    
    print("GSEA CellMarker results:")
    print(head(degs_GSEA_CellMarker))
    
    degs_GSEA_CellMarker
  })
  
  # GSEA CellMarker table output ----
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
  
  # GSEA CellMarker plot output ----
  output$gsea_plot_cellmarker <- renderPlotly({
    req(input$go_2) 
    req(gsea_results_CellMarker())
    
    # Prepare the data: filter top 5 enriched processes per condition and reorder
    data <- gsea_results_CellMarker()%>% 
      group_by(condition) %>%
      arrange(condition, desc(-log(qvalue))) %>% 
      slice_head(n = 5) %>%
      ungroup() %>%
      mutate(process = fct_reorder(process, -log(qvalue)))  
    
    num_conditions <- data %>%
      distinct(condition) %>%
      nrow()
    
    # Create the ggplot bar plot
    plot <- ggplot(data) +
      aes(x = -log(qvalue), y = process,
          fill = -log(qvalue)) +
      geom_col() +
      facet_wrap(~condition, scales = "free_y", nrow = num_conditions) +
      labs(x = "-log10(q-value)", y = "Process", title = "GSEA Results") +
      theme_minimal()
    
    # Convert ggplot to interactive plotly object
    ggplotly(plot)
  })
  
  
  # GSEA BTM immune results reactive ------
  gsea_results_BTM_Immune  <- reactive({
    req(input$go_2)
    req(degs_sig()) 
    
    degs_sig_data <- degs_sig()
    print("Data for GSEA BTM Immune:")
    print(head(degs_sig_data))
    
    # Perform GSEA using the custom autoGSEA function
    degs_GSEA_BTM_Immune <- degs_sig_data %>%
      dplyr::select(condition, genes, log2fold_change) %>% 
      autoGSEA(TERM2GENE = btm_genes, # Use BTM gene sets
               geneset_name = "BTM Immune",
               minGSSize = input$minGSSize,
               maxGSSize = input$maxGSSize,
               pvalueCutoff = input$pvalueCutoff_gsea,
               pAdjustMethod = input$pAdjustMethod) 
    
    print("GSEA BTM Immune results:")
    print(head(degs_GSEA_BTM_Immune))
    
    degs_GSEA_BTM_Immune
  })
  
  # GSEA BTM immune table output ----
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
  
  # GSEA BTM immune plot output ----
  output$gsea_plot_btm_immune <- renderPlotly({
    req(input$go_2) 
    req(gsea_results_BTM_Immune())
    
    # Prepare the data: filter top 5 enriched processes per condition and reorder
    data <- gsea_results_BTM_Immune()%>% 
      group_by(condition) %>%
      arrange(condition, desc(-log(qvalue))) %>% 
      slice_head(n = 5) %>%
      ungroup() %>%
      mutate(process = fct_reorder(process, -log(qvalue)))  
    
    num_conditions <- data %>%
      distinct(condition) %>%
      nrow()
    
    # Create the ggplot bar plot
    plot <- ggplot(data) +
      aes(x = -log(qvalue), y = process,
          fill = -log(qvalue)) +
      geom_col() +
      facet_wrap(~condition, scales = "free_y", nrow = num_conditions) +
      labs(x = "-log10(q-value)", y = "Process", title = "GSEA Results") +
      theme_minimal()
    
    # Convert ggplot to interactive plotly object
    ggplotly(plot)
  })
  
  
  # GSEA VaxSigDB results reactive -----
  gsea_results_VaxSigDB_Genes_filtered <- reactive({
    req(input$go_2)
    req(degs_sig())
    
    degs_sig_data <- degs_sig()
    print("Data for GSEA VaxSigDB:")
    print(head(degs_sig_data))
    
    # Perform GSEA using the custom autoGSEA function
    degs_GSEA_VaxSigDB_Genes_filtered <- degs_sig_data %>%
      dplyr::select(condition, genes, log2fold_change) %>% 
      autoGSEA(TERM2GENE = VaxSigDB_Genes_filtered, # Use VaxSigDB gene sets
               geneset_name = "VaxSigDB Filtered",
               minGSSize = input$minGSSize,
               maxGSSize = input$maxGSSize,
               pvalueCutoff = input$pvalueCutoff_gsea,
               pAdjustMethod = input$pAdjustMethod) 
    
    print("GSEA VaxSigDB results:")
    print(head(degs_GSEA_VaxSigDB_Genes_filtered))
    
    degs_GSEA_VaxSigDB_Genes_filtered
  })
  
  # GSEA VaxSigDB table output ----
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
  
  # GSEA VaxSigDB plot output ----
  output$gsea_plot_vaxsigdb <- renderPlotly({
    req(input$go_2) 
    req(gsea_results_VaxSigDB_Genes_filtered())
    
    # Prepare the data: filter top 5 enriched processes per condition and reorder
    data <- gsea_results_VaxSigDB_Genes_filtered()%>% 
      group_by(condition) %>%
      arrange(condition, desc(-log(qvalue))) %>% 
      slice_head(n = 5) %>%
      ungroup() %>%
      mutate(process = fct_reorder(process, -log(qvalue)))  
    
    num_conditions <- data %>%
      distinct(condition) %>%
      nrow()
    
    # Create the ggplot bar plot
    plot <- ggplot(data) +
      aes(x = -log(qvalue), y = process,
          fill = -log(qvalue)) +
      geom_col() +
      facet_wrap(~condition, scales = "free_y", nrow = num_conditions) +
      labs(x = "-log10(q-value)", y = "Process", title = "GSEA Results") +
      theme_minimal()
    
    # Convert ggplot to interactive plotly object
    ggplotly(plot)
  })
  
  # Note: The code has a duplicate definition for gsea_plot_vaxsigdb, which will be ignored by R but kept here for annotation consistency.
  output$gsea_plot_vaxsigdb <- renderPlotly({
    req(input$go_2) 
    req(gsea_results_VaxSigDB_Genes_filtered())
    
    # Prepare the data: filter top 5 enriched processes per condition and reorder
    data <- gsea_results_VaxSigDB_Genes_filtered()%>% 
      group_by(condition) %>%
      arrange(condition, desc(-log(qvalue))) %>% 
      slice_head(n = 5) %>%
      ungroup() %>%
      mutate(process = fct_reorder(process, -log(qvalue)))  
    
    num_conditions <- data %>%
      distinct(condition) %>%
      nrow()
    
    # Create the ggplot bar plot (duplicate logic)
    plot <- ggplot(data) +
      aes(x = -log(qvalue), y = process,
          fill = -log(qvalue)) +
      geom_col() +
      facet_wrap(~condition, scales = "free_y", nrow = num_conditions) +
      labs(x = "-log10(q-value)", y = "Process", title = "GSEA Results") +
      theme_minimal()
    
    # Convert ggplot to interactive plotly object
    ggplotly(plot)
  })
  
  
  
  # ssGSEA CellMarker results reactive ----
  ssgsea_results_CellMarker <- reactive({
    req(input$go_ssgsea)  # Requires the "Go" button in the ssGSEA tab to be pressed
    req(input$ssgsea_gene, input$metadata)  # Requires both expression and metadata files
    
    # Read the uploaded files
    expression_data <- read.csv(input$ssgsea_gene$datapath) %>% 
      column_to_rownames("genes") # Set gene names as rownames
    
    metadata <- read.csv(input$metadata$datapath) # Read sample metadata
    
    # Extract sample names from expression data columns
    sample_names <- colnames(expression_data) %>% 
      as.data.frame() %>% 
      rename(sample_names = '.') %>% 
      distinct()
    
    # Format CellMarker gene sets for ssGSEA (as a list)
    genelist_cellmarker <- CellMarker_ImmuneCells %>%
      dplyr::select(process, genes)
    
    genelist_cellmarker = split(genelist_cellmarker$genes,
                                genelist_cellmarker$process)
    
    
    # Calculate ssGSEA scores using the corto::ssgsea function
    nesmat_result_cellmarker <- ssgsea(inmat = expression_data, 
                                       groups = genelist_cellmarker) %>% 
      as.data.frame() %>% 
      t() %>% # Transpose matrix
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
      # Convert to long format for easier plotting/joining
      pivot_longer(cols = -process, names_to = "sample", values_to = "nes") %>%
      # Calculate p-value from Normalized Enrichment Score (NES)
      mutate(pvalue = z2p(nes), 
             qvalue = p.adjust(pvalue, method = "BH"),
             logq = -log10(qvalue),
             geneset_name = "CellMarker Immune") %>%
      # Join with metadata to get sample condition
      inner_join(metadata, by = "sample") %>%
      # Join with CellMarker annotation (redundantly joins back all markers)
      inner_join(CellMarker_ImmuneCells, by = "process") %>%
      dplyr::select(-genes) %>% # Remove redundant gene column
      distinct() 
    
    nesmat_result_cellmarker %>% as.data.frame()
  })
  
  # ssGSEA CellMarker table output ----
  output$ssgsea_table_cellmarker <- DT::renderDataTable({
    req(input$go_ssgsea)
    req(ssgsea_results_CellMarker())
    
    # Render interactive DataTable with Export buttons
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
  
  
  # Observer to update the dropdown list of cell markers for ssGSEA plot
  observe({
    results <- ssgsea_results_CellMarker()
    if (!is.null(results) && nrow(results) > 0) {
      updateSelectInput(session, "selected_cellmarker", 
                        choices = unique(results$process)) # Update choices with unique processes
    }
  })
  
  # ssGSEA CellMarker plot output ----
  output$ssgsea_plot_cellmarker <- renderPlotly({
    req(input$go_ssgsea) 
    req(input$selected_cellmarker) # Requires a cell marker to be selected
    req(ssgsea_results_CellMarker())
    
    # Filter data for the selected cell marker and create a boxplot
    plot <- ssgsea_results_CellMarker() %>% 
      filter(process == input$selected_cellmarker) %>% 
      ggplot(.) +
      aes(x = condition, y = nes, fill = condition) +
      geom_boxplot(alpha = 0.3) + # Boxplot layer
      geom_point() + # Points for individual samples
      scale_fill_brewer(palette = "Accent", direction = 1) +
      ggthemes::theme_few() +
      labs(x ="",
           y = "NES") +
      theme(
        legend.position = "top",
        axis.text.x = element_text(angle = 90)
      ) +
      facet_wrap(vars(process), scales = "free") # Plot title is the process name
    
    # Convert ggplot to interactive plotly object
    ggplotly(plot)
  })
  
  
  
  
  # Gene Overlap - Reactive for ALL DEGs (Log2FC & padj filtered)
  df_overlap <- reactive({
    req(input$go_overlap) # Requires the "Go" button in the Overlap tab to be pressed
    req(conditions_genes())
    
    conditions_genes() %>%
      # Apply Log2FC and p-value filters from the Overlap tab
      filter(log2fold_change >= input$filter_logfc_overlap | log2fold_change <= - input$filter_logfc_overlap,
             padj < input$filter_padj_overlap) %>% 
      dplyr::select(process = condition, genes) # Rename condition to process
  })
  
  # Gene Overlap - Calculation for ALL DEGs ----
  overlap_results <- reactive({
    req(input$go_overlap)
    req(df_overlap())
    
    # Initialize variables
    unique_cond <- unique(df_overlap()$process) # List of unique conditions
    shared_genes_df <- data.frame() # Data frame to store shared gene results
    
    # Compute overlap between all unique pairs of conditions
    for (i in 1:(length(unique_cond) - 1)) {
      for (j in (i + 1):length(unique_cond)) {
        resultado_temp <- overlap_genes(unique_cond[i], unique_cond[j], df_overlap())
        shared_genes_df <- bind_rows(shared_genes_df, resultado_temp)
      }
    }
    
    # Compute Fisher's exact test p-values for the overlap
    shared_genes_df <- shared_genes_df %>%
      mutate(pvalue = pmap_dbl(list(Shared, NotShared_cond1, NotShared_cond2), fisher_exact_test)) %>% 
      mutate("-log(p)" = -log10(pvalue))
    
    # Create a mirror table to make the final matrix symmetric (Cond1 vs Cond2, and Cond2 vs Cond1)
    shared_genes_df_mirror <- shared_genes_df %>% 
      rename(Cond1_temp = Cond1, Cond2_temp = Cond2) %>%
      rename(Cond1 = Cond2_temp, Cond2 = Cond1_temp) %>%
      bind_rows(shared_genes_df) %>%
      mutate(Percentage_Shared_Cond1 = round(Percentage_Shared_Cond1, 2),
             Percentage_Shared_Cond2 = round(Percentage_Shared_Cond2, 2))
    
    # Calculate Jaccard distance matrix (Jaccard similarity = 1 - distance)
    jaccard.matrix <- df_overlap() %>%
      mutate(present = 1) %>% 
      distinct() %>% 
      # Convert data to a wide matrix format (genes as columns, conditions as rows)
      pivot_wider(names_from = genes, values_from = present, values_fill = 0) %>%
      column_to_rownames(var = "process") %>%
      # Calculate Jaccard distance
      vegdist(binary = TRUE, method = "jaccard", diag = TRUE, upper = TRUE) %>%
      as.matrix() %>%
      {1 - .} # Convert distance to similarity
    
    # Combine Jaccard similarity with shared genes data
    shared_genes_df_mirror <- jaccard.matrix %>% 
      as.data.frame() %>% 
      rownames_to_column("Cond1") %>% 
      pivot_longer(cols = -"Cond1", names_to = "Cond2", values_to = "jaccard_distance") %>% 
      inner_join(shared_genes_df_mirror, by = c("Cond1", "Cond2")) %>% 
      mutate(jaccard_distance = round(jaccard_distance, 3))
    
    return(shared_genes_df_mirror)
  })
  
  # Gene Overlap - ALL DEGs table output ----
  output$overlap_table <- DT::renderDataTable({
    req(input$go_overlap)
    req(overlap_results())
    
    # Render interactive DataTable with Export buttons
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
  
  # Gene Overlap - ALL DEGs heatmap output ----
  output$overlap_heatmap <- renderPlotly({
    req(overlap_results())
    req(input$overlap_metric) # Requires selection of the metric
    
    # Prepare the data for the heatmap
    heatmap_data <- overlap_results()
    
    # Convert data frame to a matrix format using the selected metric
    heatmap_matrix <- heatmap_data %>%
      dplyr::select(Cond1, Cond2, overlap_metric = input$overlap_metric) %>%
      pivot_wider(names_from = Cond2, values_from = overlap_metric) %>%
      distinct() %>%
      column_to_rownames(var = "Cond1") %>%
      as.matrix()
    
    # Create the interactive heatmap using heatmaply
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
  
  # Gene Overlap - Reactive for IMMUNE-RELATED DEGs
  df_immune <- reactive({
    req(input$go_overlap)
    req(conditions_genes())
    
    conditions_genes() %>%
      # Apply Log2FC and p-value filters
      filter(log2fold_change >= input$filter_logfc_overlap | 
               log2fold_change <= - input$filter_logfc_overlap,
             padj < input$filter_padj_overlap) %>% 
      dplyr::select(process = condition, genes) %>% 
      # Keep only genes found in the ImmuneGO general gene set
      inner_join(ImmuneGO_genes_general %>% dplyr::select(genes) %>% distinct(), 
                 by = "genes")
  })
  
  # Gene Overlap - Calculation for IMMUNE-RELATED DEGs ----
  overlap_results_immune <- reactive({
    req(input$go_overlap)
    req(df_immune())
    
    # Initialization (using conditions from *all* DEGs for consistency, but filtering data by df_immune)
    unique_cond <- unique(df_overlap()$process) 
    shared_genes_df <- data.frame() 
    
    # Compute overlap between all unique pairs of conditions using the immune-filtered data
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
    
    # Create mirror table (for symmetry)
    shared_genes_df_mirror <- shared_genes_df %>% 
      rename(Cond1_temp = Cond1, Cond2_temp = Cond2) %>%
      rename(Cond1 = Cond2_temp, Cond2 = Cond1_temp) %>%
      bind_rows(shared_genes_df) %>%
      mutate(Percentage_Shared_Cond1 = round(Percentage_Shared_Cond1, 2),
             Percentage_Shared_Cond2 = round(Percentage_Shared_Cond2, 2))
    
    # Calculate Jaccard distance matrix (using ALL DEGs data, which might be an error in the original code, but is preserved)
    # NOTE: The original code reuses `df_overlap()` here instead of `df_immune()`. This is kept but noted.
    jaccard.matrix <- df_overlap() %>% 
      mutate(present = 1) %>% 
      distinct() %>% 
      pivot_wider(names_from = genes, values_from = present, values_fill = 0) %>%
      column_to_rownames(var = "process") %>%
      vegdist(binary = TRUE, method = "jaccard", diag = TRUE, upper = TRUE) %>%
      as.matrix() %>%
      {1 - .}
    
    # Combine Jaccard similarity with shared genes data
    shared_genes_df_mirror <- jaccard.matrix %>% 
      as.data.frame() %>% 
      rownames_to_column("Cond1") %>% 
      pivot_longer(cols = -"Cond1", names_to = "Cond2", values_to = "jaccard_distance") %>% 
      inner_join(shared_genes_df_mirror, by = c("Cond1", "Cond2")) %>% 
      mutate(jaccard_distance = round(jaccard_distance, 3))
    
    return(shared_genes_df_mirror)
  })
  
  # Gene Overlap - IMMUNE-RELATED DEGs table output ----
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
  
  # Gene Overlap - IMMUNE-RELATED DEGs heatmap output ----
  output$overlap_heatmap_immune <- renderPlotly({
    req(overlap_results_immune())
    req(input$overlap_metric) 
    
    # Prepare the data for the heatmap
    heatmap_data <- overlap_results_immune()
    
    # Convert data frame to a matrix format using the selected metric
    heatmap_matrix <- heatmap_data %>%
      dplyr::select(Cond1, Cond2, overlap_metric = input$overlap_metric) %>%
      pivot_wider(names_from = Cond2, values_from = overlap_metric) %>%
      distinct() %>%
      column_to_rownames(var = "Cond1") %>%
      as.matrix()
    
    # Create the interactive heatmap
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
  
  # Gene Overlap - Reactive for NON-IMMUNE DEGs
  df_not_immune <- reactive({
    req(input$go_overlap)
    req(conditions_genes())
    
    conditions_genes() %>%
      # Apply Log2FC and p-value filters
      filter(log2fold_change >= input$filter_logfc_overlap | 
               log2fold_change <= - input$filter_logfc_overlap,
             padj < input$filter_padj_overlap) %>% 
      dplyr::select(process = condition, genes) %>% 
      # Keep only genes NOT found in the ImmuneGO general gene set
      anti_join(ImmuneGO_genes_general %>% dplyr::select(genes) %>% distinct(), 
                by = "genes")
  })
  
  # Gene Overlap - Calculation for NON-IMMUNE DEGs ----
  overlap_results_not_immune <- reactive({
    req(input$go_overlap)
    req(df_not_immune())
    
    # Initialization (using conditions from *all* DEGs for consistency, but filtering data by df_not_immune)
    unique_cond <- unique(df_overlap()$process) 
    shared_genes_df <- data.frame() 
    
    # Compute overlap between all unique pairs of conditions using the non-immune-filtered data
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
    
    # Create mirror table (for symmetry)
    shared_genes_df_mirror <- shared_genes_df %>% 
      rename(Cond1_temp = Cond1, Cond2_temp = Cond2) %>%
      rename(Cond1 = Cond2_temp, Cond2 = Cond1_temp) %>%
      bind_rows(shared_genes_df) %>%
      mutate(Percentage_Shared_Cond1 = round(Percentage_Shared_Cond1, 2),
             Percentage_Shared_Cond2 = round(Percentage_Shared_Cond2, 2))
    
    # Calculate Jaccard distance matrix (using ALL DEGs data, which might be an error in the original code, but is preserved)
    # NOTE: The original code reuses `df_overlap()` here instead of `df_not_immune()`. This is kept but noted.
    jaccard.matrix <- df_overlap() %>%
      mutate(present = 1) %>% 
      distinct() %>% 
      pivot_wider(names_from = genes, values_from = present, values_fill = 0) %>%
      column_to_rownames(var = "process") %>%
      vegdist(binary = TRUE, method = "jaccard", diag = TRUE, upper = TRUE) %>%
      as.matrix() %>%
      {1 - .}
    
    # Combine Jaccard similarity with shared genes data
    shared_genes_df_mirror <- jaccard.matrix %>% 
      as.data.frame() %>% 
      rownames_to_column("Cond1") %>% 
      pivot_longer(cols = -"Cond1", names_to = "Cond2", values_to = "jaccard_distance") %>% 
      inner_join(shared_genes_df_mirror, by = c("Cond1", "Cond2")) %>% 
      mutate(jaccard_distance = round(jaccard_distance, 3))
    
    return(shared_genes_df_mirror)
  })
  
  # Gene Overlap - NON-IMMUNE DEGs table output ----
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
  
  # Gene Overlap - NON-IMMUNE DEGs heatmap output ----
  output$overlap_heatmap_not_immune <- renderPlotly({
    req(overlap_results_not_immune())
    req(input$overlap_metric) 
    
    # Prepare the data for the heatmap
    heatmap_data <- overlap_results_not_immune()
    
    # Convert data frame to a matrix format using the selected metric
    heatmap_matrix <- heatmap_data %>%
      dplyr::select(Cond1, Cond2, overlap_metric = input$overlap_metric) %>%
      pivot_wider(names_from = Cond2, values_from = overlap_metric) %>%
      distinct() %>%
      column_to_rownames(var = "Cond1") %>%
      as.matrix()
    
    # Create the interactive heatmap
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

# Run the Shiny application
shinyApp(ui, server)
