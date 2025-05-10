library(shiny)
library(readr)
library(readxl)
library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)
library(plotly)
library(tidyverse)
library(gprofiler2)
library(corrplot)
library(pheatmap)
library(reticulate)
library(xtable)
library(ggeasy)
library(tools)
library(zip)
library(DT)
library(rsvg)
library(png)
library(protr)
library(r3dmol)
library(UniprotR)
library(protti)
library(KSEAapp)
library(shinycssloaders)
options(warn = -1)
workmode = "funny"

source("functions.R")

#### UI functions ####
ui <- fluidPage(
  title="Proteomics Copilot",
  
  tags$head(
    tags$link(rel = "shortcut icon", href = "icon/receptor_icon.ico")
  ),
  
  tags$head(
    tags$style(HTML("
      /* Tab background */
      .nav-tabs {
        background-color: #337ab7;
        border-bottom: 1px solid #2a5d8f;
      }
      /* Inactive tabs */
      .nav-tabs > li > a {
        color: white;
        background-color: #337ab7;
        border: 1px solid #2a5d8f;
        border-bottom-color: transparent;
      }
      /* Hover state */
      .nav-tabs > li > a:hover {
        background-color: #2a5d8f;
      }
      /* Active tab */
      .nav-tabs > li.active > a,
      .nav-tabs > li.active > a:focus {
        color: white;
        background-color: #f39c12;
        border: 1px solid #d78e10;
        border-bottom-color: #f39c12;
      }
    "))
  ),
  
  titlePanel(div(
    "Proteomics Copilot",
    style = "
        background-color: #337ab7;
        color: white;
        padding: 10px;
        border-radius: 4px;
      "
  )),
  fluidRow(
    column(width = 12,
           tabsetPanel(
#### SuperTab01 - Data ####
tabPanel("Data",
          tabsetPanel(
#### Tab01 - Data Upload ####
              tabPanel(
                "Data Upload",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    h4("Upload Files"),
                    fileInput("file", "Protein Level Data", 
                              accept = c(".csv", ".tsv", ".txt", ".xlsx")),
                    fileInput("file3", "Phospho Data", 
                              accept = c(".csv", ".tsv", ".txt", ".xlsx")),
                    fileInput("file2", "Full Report", 
                              accept = c(".csv", ".tsv", ".txt", ".xlsx")),
                    hr(),
                    h4("Collapse Options"),
                    radioButtons("collapse_option", "Select Collapse Option:",
                                 choices = c("Collapsed" = "collapsed", 
                                             "Not Collapsed" = "not_collapsed"),
                                 selected = "collapsed"),
                    numericInput("collapse_cutoff", "Collapse Cutoff:", value = 0)
                  ),
                  mainPanel(
                    h4("Protein Level Data Preview"),
                    DT::dataTableOutput("table"),
                    hr(),
                    h4("Phospho Data Preview"),
                    DT::dataTableOutput("table3"),
                    hr(),
                    h4("Full Report Preview"),
                    DT::dataTableOutput("table2")
                  )
                )
              ),
#### Tab02 - Meta annotation #####
            tabPanel(
              "Data Annotation",
              sidebarLayout(
                sidebarPanel(
                  width = 3,
                  
                  h4("Normal Data Annotation"),
                  fileInput("upload_meta", "Upload Metadata (Protein Group)", accept = c(".csv", ".xlsx", ".txt", ".tsv")),
                  h5("Is the data log2 transformed?"),
                  fluidRow(
                    column(6, actionButton("log2_yes", "Yes")),
                    column(6, actionButton("log2_no", "No"))
                  ),
                  hr(),
                  numericInput("filter_num", "Filter: At least", value = 3, min = 1),
                  selectInput("filterop1", "Value(s)", choices = c("per group", "in at least one group")),
                  actionButton("apply_filter", "Apply Filter"),
                  
                  hr(),
                  
                  h4("Phospho Data Annotation"),
                  fileInput("upload_meta2", "Upload Metadata (Phospho)", accept = c(".csv", ".xlsx", ".txt", ".tsv")),
                  h5("Is the data log2 transformed?"),
                  fluidRow(
                    column(6, actionButton("log2_yes2", "Yes")),
                    column(6, actionButton("log2_no2", "No"))
                  ),
                  hr(),
                  numericInput("filter_num2", "Filter: At least", value = 3, min = 1),
                  selectInput("filterop2", "Value(s)", choices = c("per group", "in at least one group")),
                  actionButton("apply_filter2", "Apply Filter"),
                  
                  hr(),
                  
                  h4("Color Scheme"),
                  selectInput(
                    inputId = "color_palette",
                    label = "Choose a Color Palette:",
                    choices = c("Default", "Default16", "Warm/Cold", "Black/Grey", "Yue7"),
                    selected = "Default"
                  ),
                  actionButton("reloadButton", "Reload All Plots")
                ),
                
                mainPanel(
                  h4("Normal Data Condition Setup"),
                  numericInput("num_conditions", "Number of Conditions:", value = 1, min = 1),
                  uiOutput("condition_inputs"),
                  hr(),
                  
                  h4("Annotated Normal Data"),
                  DT::dataTableOutput("displayed_data"),
                  hr(),
                  
                  h4("Phospho Data Condition Setup"),
                  numericInput("num_conditions2", "Number of Conditions:", value = 1, min = 1),
                  uiOutput("condition_inputs2"),
                  hr(),
                  
                  h4("Annotated Phospho Data"),
                  DT::dataTableOutput("displayed_data2")
                )
              )
            ),
#### Tabn1 - Data imputation ####
            tabPanel(
              "Impute Data",
              sidebarLayout(
                sidebarPanel(
                  width = 3,
                  actionButton("ImputeEVE", "Impute Data Values"),
                  hr(),
                  h4("Imputation Settings"),
                  numericInput("qn1", "q-Value:", value = 0.01),
                  numericInput("adj_stdn1", "Adjust Standard Deviation:", value = 1),
                  numericInput("seedn1", "Random Seed:", value = 1337),
                  selectInput("leveln1", 
                              "Level:", 
                              choices = c("Protein", "Phosphosite")),
                  hr(),
                  downloadButton("imputed_data_down", "Download Imputed Data")
                ),
                mainPanel(
                  h4("Imputation Plots"),
                  customSpinnerWrapper(plotOutput("Impute1"), workmode = workmode),
                  customSpinnerWrapper(plotOutput("Impute2"), workmode = workmode),
                  customSpinnerWrapper(plotOutput("Impute3"), workmode = workmode),
                  h4("Imputed Data Table"),
                  DT::dataTableOutput("imputed_data_tab")
                )
              )
            ),
#### Tab3.5 - Data Distribution ####
            tabPanel("Distribution",
            ),
)),
#### SuperTab02 - QC Pipeline #####
tabPanel("QC Pipeline",
          tabsetPanel(
#### Tab03 - Coverage Plot ####
              tabPanel(
                "Coverage Plot",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    
                    actionButton("toggle_id3", "Toggle IDs"),
                    hr(),
                    selectInput("level3", "Level:",
                                choices = c("Protein", "Peptide", "Phosphosite")),
                    hr(),
                    h4("Plot Size & Resolution"),
                    numericInput("plotWidth3", "Width (cm):", value = 20),
                    numericInput("plotHeight3", "Height (cm):", value = 10),
                    numericInput("plotDPI3", "DPI:", value = 300),
                    downloadButton("downloadCoveragePlot", "Download Plot"),
                    hr(),
                    textAreaInput("text3", "Annotation Text:", value = "",
                                  width = '100%', height = '100px'),
                    actionButton("addText3up", "Add Above"),
                    actionButton("addText3down", "Add Below")
                  ),
                  mainPanel(
                    customSpinnerWrapper(
                      plotOutput("coveragePlot", height = "600px"),
                      workmode = workmode
                    )
                  )
                )
              ),
#### Tab04 - Missing Values Plot ####
              tabPanel(
                "Missing Value Plot",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    
                    selectInput("level4", "Level:",
                                choices = c("Protein", "Peptide", "Phosphosite", "Precursor")),
                    hr(),
                    h4("Plot Size & Resolution"),
                    numericInput("plotWidth4", "Width (cm):", value = 20),
                    numericInput("plotHeight4", "Height (cm):", value = 10),
                    numericInput("plotDPI4", "DPI:", value = 300),
                    downloadButton("downloadMissValPlot", "Download Plot"),
                    hr(),
                    textAreaInput("text4", "Annotation Text:", value = "",
                                  width = '100%', height = '100px'),
                    actionButton("addText4up", "Add Above"),
                    actionButton("addText4down", "Add Below")
                  ),
                  mainPanel(
                      customSpinnerWrapper(
                        plotOutput("MissValPlot", height = "600px"),
                        workmode=workmode
                      )
                  )
                )
              ),
#### Tab05 - Histogramm ####
              tabPanel(
                "Histogram Intensity",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    
                    selectInput("level5", "Level:",
                                choices = c("Protein", "Phosphosite")),
                    hr(),
                    h4("Plot Size & Resolution"),
                    numericInput("plotWidth5", "Width (cm):", value = 20),
                    numericInput("plotHeight5", "Height (cm):", value = 10),
                    numericInput("plotDPI5", "DPI:", value = 300),
                    downloadButton("downloadHistIntPlot", "Download Plot"),
                    hr(),
                    textAreaInput("text5", "Annotation Text:", value = "",
                                  width = '100%', height = '100px'),
                    actionButton("addText5up", "Add Above"),
                    actionButton("addText5down", "Add Below")
                  ),
                  mainPanel(
                    customSpinnerWrapper(
                      plotOutput("HistIntPlot", height = "600px"),
                      workmode=workmode
                    )
                  ),
                )
              ),
#### Tab06 - Boxplot ####
              tabPanel(
                "Boxplot Intensity",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    
                    actionButton("toggle_outliersBX", "Toggle Outliers"),
                    actionButton("toggle_meanBX", "Mean/Single"),
                    actionButton("toggle_id6", "Toggle ID"),
                    hr(),
                    selectInput("level6", "Level:",
                                choices = c("Protein", "Phosphosite")),
                    hr(),
                    h4("Plot Size & Resolution"),
                    numericInput("plotWidth6", "Width (cm):", value = 20),
                    numericInput("plotHeight6", "Height (cm):", value = 10),
                    numericInput("plotDPI6", "DPI:", value = 300),
                    downloadButton("downloadBoxIntPlot", "Download Plot"),
                    hr(),
                    textAreaInput("text6", "Annotation Text:", value = "",
                                  width = '100%', height = '100px'),
                    actionButton("addText6up", "Add Above"),
                    actionButton("addText6down", "Add Below")
                  ),
                  mainPanel(
                    customSpinnerWrapper(
                      plotOutput("BoxIntPlot", height = "600px"),
                      workmode=workmode
                    )
                  ),
                )
              ),
#### Tab07 - COV Plot ####
              tabPanel(
                "Coefficient of Variation Plot",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    
                    actionButton("toggle_outliersCOV", "Toggle Outliers"),
                    hr(),
                    selectInput("level7", "Level:",
                                choices = c("Protein", "Phosphosite")),
                    hr(),
                    h4("Plot Size & Resolution"),
                    numericInput("plotWidth7", "Width (cm):", value = 20),
                    numericInput("plotHeight7", "Height (cm):", value = 10),
                    numericInput("plotDPI7", "DPI:", value = 300),
                    downloadButton("downloadCovPlot", "Download Plot"),
                    hr(),
                    textAreaInput("text7", "Annotation Text:", value = "",
                                  width = '100%', height = '100px'),
                    actionButton("addText7up", "Add Above"),
                    actionButton("addText7down", "Add Below")
                  ),
                  mainPanel(
                    customSpinnerWrapper(
                      plotOutput("CovPlot", height = "600px"),
                      workmode=workmode
                    )
                  ),
                )
              ),
#### Tab08 - PCA #### 
              tabPanel(
                "PCA Plot",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    
                    selectInput("level8", "Level:",
                                choices = c("Protein", "Phosphosite")),
                    selectInput("style8", "Select Plot Type:",
                                choices = c("PCA", "tSNE", "UMAP")),
                    hr(),
                    h4("Plot Size & Resolution"),
                    numericInput("plotWidth8", "Width (cm):", value = 20),
                    numericInput("plotHeight8", "Height (cm):", value = 10),
                    numericInput("plotDPI8", "DPI:", value = 300),
                    downloadButton("downloadPCAPlot", "Download Plot"),
                    hr(),
                    textAreaInput("text8", "Annotation Text:", value = "",
                                  width = '100%', height = '100px'),
                    actionButton("addText8up", "Add Above"),
                    actionButton("addText8down", "Add Below")
                  ),
                  mainPanel(
                    customSpinnerWrapper(
                      plotOutput("PCAPlot", height = "600px"),
                      workmode=workmode
                    )
                  ),
                )
              ),
#### Tab09 - Abundance Plot ####
              tabPanel(
                "Abundance Plot",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    selectInput("level9", "Level:", choices = c("Protein", "Phosphosite")),
                    hr(),
                    h4("Plot Size & Download"),
                    numericInput("plotWidth9", "Width (cm):", value = 20),
                    numericInput("plotHeight9", "Height (cm):", value = 10),
                    numericInput("plotDPI9",    "DPI:",         value = 300),
                    downloadButton("downloadAbPlot", "Download Plot"),
                    hr(),
                    textAreaInput("text9", "Enter text:", value = "",
                                  width = '100%', height = '150px'),
                    actionButton("addText9up",   "Add Text Above Plot"),
                    actionButton("addText9down", "Add Text Below Plot")
                  ),
                  mainPanel(
                    fluidRow(
                      column(
                        width = 6,
                        selectInput("condition9", "Choose Condition:",
                                    choices = c("All Conditions"),
                                    selected = "All Conditions")
                      ),
                      column(
                        width = 6,
                        selectizeInput("protein9", "Select Proteins:",
                                       choices = NULL, multiple = TRUE,
                                       width = '100%')
                      )
                    ),
                    hr(),
                    
                    customSpinnerWrapper(
                      plotlyOutput("abundancePlot", height = "600px"),
                      workmode=workmode
                    ),
                  )
                )
              ),
#### Tab12 - Correlation Plot ####
              tabPanel(
                "Correlation Plot",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    
                    actionButton("Change12",      "Change Display"),
                    actionButton("toggle_id12",   "Toggle ID"),
                    hr(),
                    
                    selectInput(
                      "level12",
                      "Level:",
                      choices = c("Protein", "Phosphosite")
                    ),
                    hr(),
                    
                    h4("Plot Size & Resolution"),
                    numericInput(
                      "plotWidth12",
                      "Width (cm):",
                      value = 20
                    ),
                    numericInput(
                      "plotHeight12",
                      "Height (cm):",
                      value = 10
                    ),
                    numericInput(
                      "plotDPI12",
                      "DPI:",
                      value = 300
                    ),
                    downloadButton("downloadCorrPlot", "Download Plot"),
                    hr(),
                    
                    textAreaInput(
                      "text12",
                      "Enter text:",
                      value  = "",
                      width  = '100%',
                      height = '200px'
                    ),
                    actionButton("addText12up",   "Add Text Above Plot"),
                    actionButton("addText12down", "Add Text Below Plot")
                  ),
                  
                  mainPanel(
                    h4("Correlation Plot (Pearson)"),
                    customSpinnerWrapper(
                      plotOutput("CorrPlot", height = "600px"),
                      workmode=workmode
                    )
                  )
                )
              ),
#### Tab13 - Heatmap ####
              tabPanel(
                "Heatmap",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    
                    actionButton("toggle_id13", "Toggle ID"),
                    actionButton("show13",    "Show Missing Values"),
                    actionButton("fix13",     "Fix Plot"),
                    hr(),
                    
                    selectInput(
                      "level13",
                      "Level:",
                      choices = c("Protein", "Phosphosite")
                    ),
                    hr(),
                    
                    h4("Plot Size & Resolution"),
                    numericInput(
                      "plotWidth13",
                      "Width (cm):",
                      value = 20
                    ),
                    numericInput(
                      "plotHeight13",
                      "Height (cm):",
                      value = 10
                    ),
                    numericInput(
                      "plotDPI13",
                      "DPI:",
                      value = 300
                    ),
                    hr(),
                    
                    downloadButton("downloadheatPlot", "Download Plot"),
                    hr(),
                    
                    textAreaInput(
                      "text13",
                      "Enter text:",
                      value  = "",
                      width  = "100%",
                      height = "200px"
                    ),
                    actionButton("addText13up",   "Add Text Above Plot"),
                    actionButton("addText13down", "Add Text Below Plot")
                  ),
                  
                  mainPanel(
                    customSpinnerWrapper(
                      plotOutput("heatPlot", height = "800px"),
                      workmode=workmode
                    )
                  )
                )
              ),
)),
#### SuberTab03 - Statistical Analysis ####
tabPanel("Statistical Analysis",
          tabsetPanel(
#### Tab10 - Volcano Plot ####
              tabPanel(
                "Volcano Plot",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    
                    h4("Thresholds"),
                    numericInput("in_pval10",   "P‑value threshold:",    value = 0.05, step = 0.01),
                    numericInput("in_log2fc10", "log₂ FC threshold:",    value = 1,    step = 0.1),
                    checkboxInput("uncorrected10", "Use uncorrected p‑values", value = FALSE),
                    hr(),
                    
                    selectInput("paired10", "Test Type:",   choices = c("Unpaired", "Paired")),
                    selectInput("level10",  "Level:",       choices = c("Protein", "Phosphosite")),
                    hr(),
                    
                    downloadButton("download10", "Download Data"),
                    hr(),
                    
                    textAreaInput("text10", "Annotation Text:", value = "",
                                  width = '100%', height = '100px'),
                    actionButton("addText10up",   "Add Above"),
                    actionButton("addText10down", "Add Below"),
                    hr(),
                    
                    actionButton("addVolc", "Add Plot to Report", 
                                 icon = icon("file-medical"), 
                                 style = "width:100%; margin-top: 20px;")
                  ),
                  
                  mainPanel(
                    fluidRow(
                      column(6, selectInput("condition1_10", "Condition 1:", choices = c())),
                      column(6, selectInput("condition2_10", "Condition 2:", choices = c()))
                    ),
                    hr(),
                    
                    customSpinnerWrapper(
                      plotlyOutput("VolcPlot", height = "500px"),
                      workmode=workmode
                    ),
                    
                    hr(),
                    
                    DT::dataTableOutput("table10")
                  )
                )
              ),
#### Tab11 - GSEA ####
              tabPanel(
                "GSEA",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    
                    numericInput(
                      "top_n11",
                      "Number of top gene sets to show:",
                      value = 10,
                      min   = 1,
                      step  = 1
                    ),
                    hr(),
                    
                    h4("Term Size Filter"),
                    numericInput(
                      "filter11min",
                      "Min term size:",
                      value = 20
                    ),
                    numericInput(
                      "filter11max",
                      "Max term size:",
                      value = 300
                    ),
                    hr(),
                    
                    h4("Plot Size & Resolution"),
                    numericInput(
                      "plotWidth11",
                      "Width (cm):",
                      value = 20
                    ),
                    numericInput(
                      "plotHeight11",
                      "Height (cm):",
                      value = 10
                    ),
                    numericInput(
                      "plotDPI11",
                      "DPI:",
                      value = 300
                    ),
                    hr(),
                    
                    downloadButton("downloadUpPlot",   "Download Plot (UP)"),
                    downloadButton("downloadDownPlot", "Download Plot (DOWN)")
                  ),
                  
                  mainPanel(
                    h3("Upregulated Gene Sets"),
                    customSpinnerWrapper(
                      plotOutput("UpregEnrichmentPlot", height = "600px"),
                      workmode=workmode
                    ),
                    hr(),
                    
                    h3("Downregulated Gene Sets"),
                    customSpinnerWrapper(
                      plotOutput("DownregEnrichmentPlot", height = "600px"),
                      workmode=workmode
                    ),
                    hr(),
                    
                    fluidRow(
                      column(
                        6,
                        h4("Upregulated Gene List"),
                        DT::dataTableOutput("table11up")
                      ),
                      column(
                        6,
                        h4("Downregulated Gene List"),
                        DT::dataTableOutput("table11down")
                      )
                    )
                  )
                )
              ),
)),
#### SuperTab04 - Peptide Level #####
tabPanel("Peptide Level",
         tabsetPanel(
#### Tab14 - RT Plot ####
              tabPanel(
                "RT Plot",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    
                    actionButton("line14", "Add line"),
                    hr(),
                    
                    selectInput(
                      "style14",
                      "Select Plot Type:",
                      choices = c("Scatter Plot", "Hexbin Plot", "Density Plot")
                    ),
                    numericInput(
                      "bins14",
                      "Bins:",
                      value = 1000
                    ),
                    hr(),
                    
                    h4("Plot Size & Resolution"),
                    numericInput(
                      "plotWidth14",
                      "Width (cm):",
                      value = 20
                    ),
                    numericInput(
                      "plotHeight14",
                      "Height (cm):",
                      value = 10
                    ),
                    numericInput(
                      "plotDPI14",
                      "DPI:",
                      value = 300
                    ),
                    downloadButton("downloadRTPlot", "Download Plot"),
                    hr(),
                    
                    textAreaInput(
                      "text14",
                      "Enter text:",
                      value  = "",
                      width  = '100%',
                      height = '200px'
                    ),
                    actionButton("addText14up",   "Add Text Above Plot"),
                    actionButton("addText14down", "Add Text Below Plot")
                  ),
                  
                  mainPanel(
                    customSpinnerWrapper(
                      plotOutput("RTPlot", height = "600px"),
                      workmode=workmode
                    )
                  )
                )
              ),
#### Tab15 - Modification Plot ####
              tabPanel(
                "Modification plot",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    
                    actionButton("toggle_id15", "Toggle ID"),
                    hr(),
                    
                    h4("Plot Size & Resolution"),
                    numericInput("plotWidth15",  "Width (cm):",  value = 20),
                    numericInput("plotHeight15", "Height (cm):", value = 10),
                    numericInput("plotDPI15",    "DPI:",        value = 300),
                    downloadButton("downloadModPlot", "Download Plot"),
                    hr(),
                    
                    textAreaInput(
                      "text15",
                      "Enter text:",
                      value  = "",
                      width  = "100%",
                      height = "200px"
                    ),
                    actionButton("addText15up",   "Add Text Above Plot"),
                    actionButton("addText15down", "Add Text Below Plot")
                  ),
                  
                  mainPanel(
                    customSpinnerWrapper(
                      plotOutput("ModPlot", height = "600px"),
                      workmode=workmode
                    )
                  )
                )
              ),
#### Tab16 - Missed Cleavage Plot ####
              tabPanel(
                "Missed Cleavage Plot",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    
                    actionButton("toggle_id16", "Toggle ID"),
                    hr(),
                    
                    h4("Plot Size & Resolution"),
                    numericInput("plotWidth16",  "Width (cm):",  value = 20),
                    numericInput("plotHeight16", "Height (cm):", value = 10),
                    numericInput("plotDPI16",    "DPI:",         value = 300),
                    downloadButton("downloadMCPlot", "Download Plot"),
                    hr(),
                    
                    textAreaInput(
                      "text16",
                      "Enter text:",
                      value  = "",
                      width  = "100%",
                      height = "200px"
                    ),
                    actionButton("addText16up",   "Add Text Above Plot"),
                    actionButton("addText16down", "Add Text Below Plot")
                  ),
                  
                  mainPanel(
                    customSpinnerWrapper(
                      plotOutput("MCPlot", height = "600px"),
                      workmode=workmode
                    )
                  )
                )
              ),
#### Tab22 - Protein centric ####
              tabPanel("Protein Centric",
                       actionButton("db_fasta22", "Load Database"),
                       selectizeInput("protein20", "Select Protein:", choices = NULL),
                       numericInput("chunk20", "Chunk size:", value = 150),
                       verbatimTextOutput("seq_allign"),
                       r3dmolOutput("plot_3d22"),
              ),
)),
#### SuberTab05 - Single Protein ####
tabPanel("Single Protein",
         tabsetPanel(
#### Tab17 - Boxplot, single protein ####
             tabPanel("Protein Box",
                      selectizeInput("protein17", "Select Protein:", choices = NULL),
                      selectizeInput("conditions17", "Select Conditions:", choices = NULL, multiple = TRUE),
                      plotOutput("ProtBoxPlot"),
                      selectInput("level17", 
                                  "Level:", 
                                  choices = c("Protein", "Phosphosite"))
             ),
             
#### Tab18 - Lineplot, single protein ####
             tabPanel("Protein Line",
                      selectizeInput("protein18", "Select Proteins:", choices = NULL, multiple = TRUE),
                      selectizeInput("conditions18", "Select Conditions:", choices = NULL, multiple = TRUE),
                      plotOutput("ProtLinePlot"),
                      selectInput("level18", 
                                  "Level:", 
                                  choices = c("Protein", "Phosphosite"))
             ),
)),
#### SuperTab06 - Phospho-specific #####
tabPanel("Phospho-specific",
         tabsetPanel(
#### Tab19 - Phossite Plot ####
            tabPanel(
              "Phossite Plot",
              sidebarLayout(
                sidebarPanel(
                  width = 3,
                  numericInput("cutoff19", "Cutoff:", value = 0),
                  hr(),
                  h4("Annotation Text"),
                  textAreaInput("text19", "Enter text:", value = "", 
                                width = '100%', height = '100px'),
                  actionButton("addText19up", "Add Text Above Plot"),
                  actionButton("addText19down", "Add Text Below Plot")
                ),
                mainPanel(
                  h4("Phossite Plot"),
                  customSpinnerWrapper(
                    plotOutput("PhossitePlot", height = "600px"),
                    workmode=workmode
                  )
                )
              )
            ),
#### Tab21 - Phossite Coverage Plot ####
          tabPanel(
            "Phossite Coverage Plot",
            sidebarLayout(
              sidebarPanel(
                width = 3,
                actionButton("toggle_id21", "Toggle IDs"),
                hr(),
                h4("Plot Size & Resolution"),
                numericInput("plotWidth21", "Width (cm):", value = 20),
                numericInput("plotHeight21", "Height (cm):", value = 10),
                numericInput("plotDPI21", "DPI:", value = 300),
                downloadButton("downloadPhossitePlot", "Download Plot"),
                hr(),
                h4("Annotation Text"),
                textAreaInput("text21", "Enter text:", value = "", width = '100%', height = '100px'),
                actionButton("addText21up", "Add Text Above Plot"),
                actionButton("addText21down", "Add Text Below Plot")
              ),
              mainPanel(
                h4("Phossite Coverage Plot"),
                customSpinnerWrapper(
                  plotOutput("PhossiteCoveragePlot", height = "600px"),
                  workmode=workmode
                )
              )
            )
          ),
#### Tab20.1 - KSEA ####
              tabPanel("KSEA",
                       fluidRow(
                         column(width = 6,
                                selectInput("condition1_20", "Select Condition 1:", 
                                            choices = c())
                         ),
                         column(width = 6,
                                selectInput("condition2_20", "Select Condition 2:", 
                                            choices = c())
                         )
                       ),
                       downloadButton("downKSEA", "Download KSEA input data"),
                       br(),
                       br(),
                       p("1. Select the conditions you want to compare"),
                       p("2. Download the data table"),
                       p("3. Visit the following website:"),
                       tags$a(href = "https://www.phosphosite.org/kinaseLibraryAction", "Website for KSEA", target = "_blank"),
                       br(),
                       br(),
                       p("4. Go to the tab >Enrichment analysis<, upload your data and click >Run enrichment analysis<"),
                       p("5. Download the table"),
                       fileInput("KSEA_data", "Upload the file here", 
                                 accept = c(".txt")),
                       plotlyOutput("KinaseVolcPlot"),
                       numericInput("in_pval20", "P-value threshold: ", value = 0.1),
                       downloadButton("downTREE", "Download Kinase Tree input data"),
                       br(),
                       br(),
                       p("6. Visit the following website:"),
                       tags$a(href = "http://phanstiel-lab.med.unc.edu/CORAL/", "Website for Kinase Trees", target = "_blank"),
                       br(),
                       br(),
                       p("7. On the website:"),
                       p("Branch Color > Color Scheme > Select Quantative > Kinases & Value > Paste enrichment_scores.txt values"),
                       p("Node Color > Color Scheme > Select Quantative > Kinases & Value > Paste enrichment_scores.txt values"),
                       p("Node Size > Scaling Scheme > Select Quantative > Kinases & Value > Paste pvals.txt values"),
                       p("Node Size > Missing Kinases > Tick hide"),
                       fileInput("TreeSVG", "Upload and Display Kinase Tree SVG Image", accept = c("image/svg+xml")),
                       uiOutput("treePlot")
              ),

#### Tab20.2 - KSEA (Kinact) ####
              tabPanel("KSEA (Kinact)",
                       fluidRow(
                         column(width = 6,
                                selectInput("condition1_202", "Select Condition 1:", 
                                            choices = c())
                         ),
                         column(width = 6,
                                selectInput("condition2_202", "Select Condition 2:", 
                                            choices = c())
                         )
                       ),
                       plotlyOutput("KinactPlot"),
                       fluidRow(
                         column(width = 6,
                                numericInput("top_n202", "Display Top or Bottom n Kinases: ", value = 0)
                         ),
                         column(width = 6,
                                numericInput("m.cutoff202", "Minimal number of Phosphosites identified in your Dataset for a Kinase: ", value = 5)
                         )
                       ),
                       fluidRow(
                         column(width = 6,
                                selectInput("NetworKIN202", 
                                            "Include NetworKIN prediction?", 
                                            choices = c(FALSE, TRUE))
                         ),
                         column(width = 6,
                                numericInput("NetworKIN.cutoff202", "Minimun NetworKIN score: ", value = 1)
                         )
                       ),
                       textInput("Kinase_202", "Enter a Kinase:"),
                       plotlyOutput("KinactPlotSites"),
                ),

)),
#### SuperTab07 - Tables ####
tabPanel("Tables",
         tabsetPanel(            
#### TabSummary - Summary ####
             tabPanel("Summary", 
                      h4("Data object"),
                      DT::dataTableOutput("transformed_data_prot"),
                      h4("Meta object"),
                      DT::dataTableOutput("meta_object"),
                      downloadButton("download_meta", "Download Metadata"),
                      h4("Phos data object"),
                      DT::dataTableOutput("collapse_data_sum"),
                      downloadButton("download_collapse_data_sum", "Download Collapse Data"),
                      h4("Meta object phos"),
                      DT::dataTableOutput("meta_object2"),
                      downloadButton("download_meta2", "Download Metadata Phos"),
             ),
             
#### Log ####
             #Log
             tabPanel("Log",
                      downloadButton("log", "Download Log"),    
                      tableOutput("LogTable"),
             ),
)),
      )
    )
  )
)

#### Server functions ####
server <- function(input, output, session) {
  options(shiny.maxRequestSize=10*1024*1024^2)
#### Tab01 - Data upload ####
  data <- reactive({
    req(input$file)
    
    ext <- tools::file_ext(input$file$name)
    df <- switch(ext,
                 csv = read_csv(input$file$datapath),
                 tsv = read_tsv(input$file$datapath),
                 txt = read_delim(input$file$datapath, delim = "\t"),
                 xlsx = read_excel(input$file$datapath),
                 stop("Invalid file type")
    )
    
    df = rename_cols(df)
    return(df)
  })
  
  data2 <- reactive({
    req(input$file2)
    
    ext <- tools::file_ext(input$file2$name)
    df <- switch(ext,
                 csv = read_csv(input$file2$datapath),
                 tsv = read_tsv(input$file2$datapath),
                 txt = read_delim(input$file2$datapath, delim = "\t"),
                 xlsx = read_excel(input$file2$datapath),
                 stop("Invalid file type")
    )
    
    df = rename_cols(df)
    return(df)
  })
  
  data3 <- reactive({
    req(input$file3)
    
    ext <- tools::file_ext(input$file3$name)
    df <- switch(ext,
                 csv = read_csv(input$file3$datapath),
                 tsv = read_tsv(input$file3$datapath),
                 txt = read_delim(input$file3$datapath, delim = "\t"),
                 xlsx = read_excel(input$file3$datapath),
                 stop("Invalid file type")
    )
    if (input$collapse_option == "not_collapsed"){
      python_script = "Collapse.py"
      command <- sprintf("python %s %s %s", python_script, input$file3$name, input$collapse_cutoff)
      system(command)
      df = read.csv(paste0(getwd(), "/Data/", file_path_sans_ext(input$file3$name), "_collapsed", ".csv"))
    }
    return(df)
  })
  
  output$table <- DT::renderDataTable({
    req(data())
    DT::datatable(data(), options = list(pageLength = 5))
  })
  
  output$table2 <- DT::renderDataTable({
    req(data2())
    DT::datatable(data2(), options = list(pageLength = 5))
  })
  
  output$table3 <- DT::renderDataTable({
    req(data3())
    DT::datatable(data3(), options = list(pageLength = 5))
  })
  
#### Tab02 - Data annotation ####
  #Normal
  output$condition_inputs <- renderUI({
    req(data())
    num_conditions <- input$num_conditions
    lapply(1:num_conditions, function(i) {
      fluidRow(
        column(width = 6,
               textInput(inputId = paste0("condition_", i), label = paste("Condition", i))
        ),
        column(width = 6,
               selectizeInput(
                 inputId = paste0("columns_", i),
                 label = paste("Select columns for Condition", i),
                 choices = colnames(data()),
                 selected = NULL,
                 multiple = TRUE
               )
        )
      )
    })
  })
  
  observe({
    lapply(1:input$num_conditions, function(i) {
      condition_id <- paste0("condition_", i)
      observeEvent(input[[condition_id]], {
        updateSelectizeInput(session, paste0("columns_", i), label = paste("Select columns for", input[[condition_id]]))
      })
    })
  })
  
  meta <- reactive({
    req(data())
    
    if (!is.null(input$upload_meta)) {
      file_ext <- tools::file_ext(input$upload_meta$name)
      
      if (file_ext == "csv") {
        meta_df <- read.csv(input$upload_meta$datapath, stringsAsFactors = FALSE)
      } else if (file_ext == "xlsx") {
        meta_df <- read_excel(input$upload_meta$datapath)
      } else if (file_ext == "txt") {
        meta_df <- read_delim(input$upload_meta$datapath, delim = "\t", col_types = cols())
      } else if (file_ext == "tsv") {
        meta_df <- read_data(input$upload_meta$datapath)
      } else {
        stop("Wrong file type!")
      }
      
      if (grepl("Setup", input$upload_meta$name)) {
        meta_df <- transform_meta(data(), meta_df)
      }
      
      return(meta_df)
    }
    
    conditions <- lapply(1:input$num_conditions, function(i) {
      list(
        condition = input[[paste0("condition_", i)]],
        columns = input[[paste0("columns_", i)]]
      )
    })
    
    meta_df <- do.call(rbind, lapply(conditions, function(cond) {
      data.frame(
        sample = unlist(cond$columns),
        condition = cond$condition,
        stringsAsFactors = FALSE
      )
    }))
    
    return(meta_df)
  })
  
  was_transformed <- reactiveVal(FALSE)
  transformed_data <- reactiveVal(NULL)
  not_transformed_data <- reactiveVal(NULL)
  
  observeEvent(input$log2_yes, {
    tryCatch({
      req(data(), meta())
      transformed_data(data())
      not_transformed_data(inverseof_log2_transform_data(data(), meta()))
      was_transformed(FALSE)
    }, error = function(e) {
      showNotification("Meta data is missing or invalid. Please define it before applying log2 settings.", type = "error")
    })
  })
  
  observeEvent(input$log2_no, {
    tryCatch({
      req(data(), meta())
      transformed_data(log2_transform_data(data(), meta()))
      not_transformed_data(data())
      was_transformed(TRUE)
    }, error = function(e) {
      showNotification("Meta data is missing or invalid. Please define it before applying log2 settings.", type = "error")
    })
  })
  
  
  filtered_data <- reactiveVal(NULL)
  log2_filtered_data <- reactiveVal(NULL)
  
  observeEvent(input$apply_filter, {
    req(data(), transformed_data(), meta(), input$filter_num)
    filtered_data(filter_data(data(), meta(), num = input$filter_num, filterops = input$filterop1))
    log2_filtered_data(filter_data(transformed_data(), meta(), num = input$filter_num, filterops = input$filterop1))
  })
  
  
  output$displayed_data <- DT::renderDataTable({
    if (is.null(log2_filtered_data())) {
      req(transformed_data())
      DT::datatable(transformed_data(), options = list(pageLength = 5))
    } else {
      req(log2_filtered_data())
      DT::datatable(log2_filtered_data(), options = list(pageLength = 5))
    }
  })
  
  volcano_data <- reactive({
    if (!is.null(imputed_data())) {
      req(imputed_data())
      return(imputed_data())
    } else {
      if (!is.null(log2_filtered_data())) {
        return(log2_filtered_data())
      } else {
        req(transformed_data())
        return(transformed_data())
      }
    }
  })
  
  #Phospho
  output$condition_inputs2 <- renderUI({
    req(data3())
    num_conditions2 <- input$num_conditions2
    lapply(1:num_conditions2, function(i) {
      fluidRow(
        column(width = 6,
               textInput(inputId = paste0("condition2_", i), label = paste("Condition", i))
        ),
        column(width = 6,
               selectizeInput(
                 inputId = paste0("columns2_", i),
                 label = paste("Select columns for Condition", i),
                 choices = colnames(data3()),
                 selected = NULL,
                 multiple = TRUE
               )
        )
      )
    })
  })
  
  observe({
    lapply(1:input$num_conditions2, function(i) {
      condition_id2 <- paste0("condition2_", i)
      observeEvent(input[[condition_id2]], {
        updateSelectizeInput(session, paste0("columns2_", i), label = paste("Select columns for", input[[condition_id2]]))
      })
    })
  })
  
  meta2 <- reactive({
    req(data3())
    
    if (!is.null(input$upload_meta2)) {
      file_ext <- tools::file_ext(input$upload_meta2$name)
      
      if (file_ext == "csv") {
        meta_df2 <- read.csv(input$upload_meta2$datapath, stringsAsFactors = FALSE)
      } else if (file_ext == "xlsx") {
        meta_df2 <- read_excel(input$upload_meta2$datapath)
      } else if (file_ext == "txt") {
        meta_df2 <- read_delim(input$upload_meta2$datapath, delim = "\t", col_types = cols())
      } else if (file_ext == "txt") {
        meta_df2 <- read_data(input$upload_meta2$datapath)
      } else {
        stop("Wrong file type!")
      }
      
      if (grepl("Setup", input$upload_meta2$name)) {
        meta_df2 <- transform_meta(data(), meta_df2)
      }
      
      return(meta_df2)
    }
    
    conditions2 <- lapply(1:input$num_conditions2, function(i) {
      list(
        condition = input[[paste0("condition2_", i)]],
        columns = input[[paste0("columns2_", i)]]
      )
    })
    
    meta_df2 <- do.call(rbind, lapply(conditions2, function(cond) {
      data.frame(
        sample = unlist(cond$columns),
        condition = cond$condition,
        stringsAsFactors = FALSE
      )
    }))
    
    return(meta_df2)
  })
  
  was_transformed3 <- reactiveVal(FALSE)
  transformed_data3 <- reactiveVal(NULL)
  not_transformed_data3 <- reactiveVal(NULL)
  
  observeEvent(input$log2_yes2, {
    tryCatch({
      req(data3(), meta2())
      transformed_data3(data3())
      not_transformed_data3(inverseof_log2_transform_data(data3(), meta2()))
      was_transformed3(FALSE)
    }, error = function(e) {
      showNotification("Meta data is missing or invalid. Please define it before applying log2 settings.", type = "error")
    })
  })
  
  observeEvent(input$log2_no2, {
    tryCatch({
      req(data3(), meta2())
      transformed_data3(log2_transform_data(data3(), meta2()))
      not_transformed_data3(data3())
      was_transformed3(TRUE)
    }, error = function(e) {
      showNotification("Meta data is missing or invalid. Please define it before applying log2 settings.", type = "error")
    })
  })
  
  filtered_data3 <- reactiveVal(NULL)
  log2_filtered_data3 <- reactiveVal(NULL)

  observeEvent(input$apply_filter2, {
    req(data3(), transformed_data3(), meta2(), input$filter_num2)
    filtered_data3(filter_data(data3(), meta2(), num = input$filter_num2, filterops = input$filterop2))
    log2_filtered_data3(filter_data(transformed_data3(), meta2(), num = input$filter_num2, filterops = input$filterop2))
  })  
  
  output$displayed_data2 <- DT::renderDataTable({
    if (is.null(log2_filtered_data3())) {
      req(transformed_data3())
      DT::datatable(transformed_data3(), options = list(pageLength = 5))
    } else {
      req(log2_filtered_data3())
      DT::datatable(log2_filtered_data3(), options = list(pageLength = 5))
    }
  })
  
  volcano_data3 <- reactive({
    if (!is.null(log2_filtered_data3())) {
      return(log2_filtered_data3())
    } else {
      req(transformed_data3())
      return(transformed_data3())
    }
  })
  
  #Color
  observeEvent(input$color_palette, {
    init_colors(color = input$color_palette)
  })
  
  observeEvent(input$reloadButton, {
    #CoveragePlot
    output$coveragePlot <- renderPlot({
      if (input$level3 == "Protein") {
        req(data(), meta())
        coverage_plot(data(), meta(), id3())
      } else if (input$level3 == "Peptide") {  
        req(data2(), meta())
        coverage_plot_pep(data2(), meta(), id3())
      } else if (input$level3 == "Phosphosite") {  
        req(data3(), meta2())
        coverage_plot(data3(), meta2(), id3())
      }
    })
    #HistoInt
    output$HistIntPlot <- renderPlot({
      if (input$level5 == "Protein") {
        req(transformed_data(), meta())
        histo_int(transformed_data(), meta())
      } else if (input$level5 == "Phosphosite"){
        req(data3(), meta2())
        histo_int(data3(), meta2())
      }
    })
    #BoxplotInt
    output$BoxIntPlot <- renderPlot({
      if (input$level6 == "Protein") {
        req(transformed_data(), meta())
        if (meanBX()){
          boxplot_int(transformed_data(), meta(), outliers = outliersBX())
        }else{
          boxplot_int_single(transformed_data(), meta(), outliers = outliersBX(), id=id6()) 
        }
      } else if (input$level6 == "Phosphosite"){
        req(data3(), meta2())
        if (meanBX()){
          boxplot_int(data3(), meta2(), outliers = outliersBX())
        }else{
          boxplot_int_single(data3(), meta2(), outliers = outliersBX(), id=id6()) 
        }
      }
    })
    #Cov Plot
    output$CovPlot <- renderPlot({
      if (input$level7 == "Protein") {
        req(data(), meta())
        cov_plot(not_transformed_data(), meta(), outliers = outliersCOV())
      } else if (input$level7 == "Phosphosite"){
        req(data3(), meta2())
        cov_plot(not_transformed_data3(), meta2(), outliers = outliersCOV())
      }
    })
    #PCA Plot
    output$PCAPlot <- renderPlot({
      if (input$level8 == "Protein") {
        req(data(), meta())
        dim_func(data(), meta(), method=input$style8)
      } else if (input$level8 == "Phosphosite"){
        req(data3(), meta2())
        dim_func(data3(), meta2(), method=input$style8)
      }
    })
    #Abundance Plot
    output$abundancePlot <- renderPlotly({
      if (input$level9 == "Protein"){
        req(input$condition9)
        if (input$condition9 == "All Conditions") {
          abundance_plot(not_transformed_data(), meta(), workflow = input$level9)
        } else {
          interactive_abundance_plot(not_transformed_data(), meta(), input$condition9, workflow = input$level9, search = input$protein9)
        }} else if (input$level9 == "Phosphosite"){
          req(input$condition9)
          if (input$condition9 == "All Conditions") {
            abundance_plot(not_transformed_data3(), meta2(), workflow = "Phosphosite")
          } else {
            interactive_abundance_plot(not_transformed_data3(), meta2(), input$condition9, workflow = "Phosphosite", search = input$protein9)
          }}
    })
  })
  
#### Tab03 - Coverage Plot ####
  id3 <- reactiveVal(TRUE)
  
  observeEvent(input$toggle_id3, {
    id3(!id3())
  })
  
  output$coveragePlot <- renderPlot({
    if (input$level3 == "Protein") {
      req(data(), meta())
      coverage_plot(data(), meta(), id3())
    } else if (input$level3 == "Peptide") {  
      req(data2(), meta())
      coverage_plot_pep(data2(), meta(), id3())
    } else if (input$level3 == "Phosphosite") {  
      req(data3(), meta2())
      coverage_plot(data3(), meta2(), id3())
    }
  })
  
  output$downloadCoveragePlot <- downloadHandler(
    filename = function() {
      paste("coverage_plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      width <- input$plotWidth3
      height <- input$plotHeight3
      dpi <- input$plotDPI3
      
      ggsave(file, plot = coverage_plot(data(), meta()), device = "png", 
             width = width, height = height, units = "cm", dpi = dpi)
    }
  )
  
  observeEvent(input$addText3up, {
    df <- const_df()
    if (input$level3 == "Protein"){
      df <- df[df$Var != "text3upProtein", ]
      new_row <- data.frame(Var = "text3upProtein", Select = input$text3, stringsAsFactors = FALSE)
    } else if (input$level3 == "Peptide"){
      df <- df[df$Var != "text3upPeptide", ]
      new_row <- data.frame(Var = "text3upPeptide", Select = input$text3, stringsAsFactors = FALSE)
    } else if (input$level3 == "Phosphosite"){
      df <- df[df$Var != "text3upPhosphosite", ]
      new_row <- data.frame(Var = "text3upPhosphosite", Select = input$text3, stringsAsFactors = FALSE)
    }
    updateTextAreaInput(session, "text3", value = "")
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$addText3down, {
    df <- const_df()
    if (input$level3 == "Protein"){
      df <- df[df$Var != "text3downProtein", ]
      new_row <- data.frame(Var = "text3downProtein", Select = input$text3, stringsAsFactors = FALSE)
    } else if (input$level3 == "Peptide"){
      df <- df[df$Var != "text3downPeptide", ]
      new_row <- data.frame(Var = "text3downPeptide", Select = input$text3, stringsAsFactors = FALSE)
    } else if (input$level3 == "Phosphosite"){
      df <- df[df$Var != "text3downPhosphosite", ]
      new_row <- data.frame(Var = "text3downPhosphosite", Select = input$text3, stringsAsFactors = FALSE)
    }
    updateTextAreaInput(session, "text3", value = "")
    const_df(rbind(df, new_row))
  })
  
#### Tab04 - Missing Value Plot ####
  output$MissValPlot <- renderPlot({
    if (input$level4 == "Protein") {
      req(data(), meta())
      missing_value_plot(data(), meta())
    } else if (input$level4 == "Precursor"){
      req(data2(), meta())
      missing_value_plot_prec(data2(), meta())
    } else if (input$level4 == "Peptide"){
      req(data2(), meta())
      missing_value_plot_pep(data2(), meta())
    } else if (input$level4 == "Phosphosite"){
      req(data3(), meta2())
      missing_value_plot(data3(), meta2())
    }
  })
  
  output$downloadMissValPlot <- downloadHandler(
    filename = function() {
      paste("missval_plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      width <- input$plotWidth4
      height <- input$plotHeight4
      dpi <- input$plotDPI4
      
      ggsave(file, plot = missing_value_plot(data(), meta()), device = "png", 
             width = width, height = height, units = "cm", dpi = dpi)
    }
  )
  
  observeEvent(input$addText4up, {
    df <- const_df()
    if (input$level4 == "Protein"){
      df <- df[df$Var != "text4upProtein", ]
      new_row <- data.frame(Var = "text4upProtein", Select = input$text4, stringsAsFactors = FALSE)
    } else if (input$level4 == "Peptide"){
      df <- df[df$Var != "text4upPeptide", ]
      new_row <- data.frame(Var = "text4upPeptide", Select = input$text4, stringsAsFactors = FALSE)
    } else if (input$level4 == "Phosphosite"){
      df <- df[df$Var != "text4upPhosphosite", ]
      new_row <- data.frame(Var = "text4upPhosphosite", Select = input$text4, stringsAsFactors = FALSE)
    }
    updateTextAreaInput(session, "text4", value = "")
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$addText4down, {
    df <- const_df()
    if (input$level4 == "Protein"){
      df <- df[df$Var != "text4downProtein", ]
      new_row <- data.frame(Var = "text4downProtein", Select = input$text4, stringsAsFactors = FALSE)
    } else if (input$level4 == "Peptide"){
      df <- df[df$Var != "text4downPeptide", ]
      new_row <- data.frame(Var = "text4downPeptide", Select = input$text4, stringsAsFactors = FALSE)
    } else if (input$level4 == "Phosphosite"){
      df <- df[df$Var != "text4downPhosphosite", ]
      new_row <- data.frame(Var = "text4downPhosphosite", Select = input$text4, stringsAsFactors = FALSE)
    }
    updateTextAreaInput(session, "text4", value = "")
    const_df(rbind(df, new_row))
  })
  
#### Tab05 - Histogram ####
  output$HistIntPlot <- renderPlot({
    if (input$level5 == "Protein") {
      req(transformed_data(), meta())
      histo_int(transformed_data(), meta())
    } else if (input$level5 == "Phosphosite"){
      req(data3(), meta2())
      histo_int(data3(), meta2())
    }
  })
  
  output$downloadHistIntPlot <- downloadHandler(
    filename = function() {
      paste("histint_plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      width <- input$plotWidth5
      height <- input$plotHeight5
      dpi <- input$plotDPI5
      
      ggsave(file, plot = histo_int(transformed_data(), meta()), device = "png", 
             width = width, height = height, units = "cm", dpi = dpi)
    }
  )
  
  observeEvent(input$addText5up, {
    df <- const_df()
    if (input$level5 == "Protein"){
      df <- df[df$Var != "text5upProtein", ]
      new_row <- data.frame(Var = "text5upProtein", Select = input$text5, stringsAsFactors = FALSE)
    } else if (input$level5 == "Phosphosite"){
      df <- df[df$Var != "text5upPhosphosite", ]
      new_row <- data.frame(Var = "text5upPhosphosite", Select = input$text5, stringsAsFactors = FALSE)
    }
    updateTextAreaInput(session, "text5", value = "")
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$addText5down, {
    df <- const_df()
    if (input$level5 == "Protein"){
      df <- df[df$Var != "text5downProtein", ]
      new_row <- data.frame(Var = "text5downProtein", Select = input$text5, stringsAsFactors = FALSE)
    } else if (input$level5 == "Phosphosite"){
      df <- df[df$Var != "text5downPhosphosite", ]
      new_row <- data.frame(Var = "text5downPhosphosite", Select = input$text5, stringsAsFactors = FALSE)
    }
    updateTextAreaInput(session, "text5", value = "")
    const_df(rbind(df, new_row))
  })
  
#### Tab06 - Boxplot ####
  id6 <- reactiveVal(TRUE)
  
  observeEvent(input$toggle_id6, {
    id6(!id6())
  })
  
  outliersBX <- reactiveVal(FALSE)
  meanBX <- reactiveVal(FALSE)
  
  observeEvent(input$toggle_outliersBX, {
    outliersBX(!outliersBX())
  })
  
  observeEvent(input$toggle_meanBX, {
    meanBX(!meanBX())
  })
  
  output$BoxIntPlot <- renderPlot({
    if (input$level6 == "Protein") {
      req(transformed_data(), meta())
      if (meanBX()){
        boxplot_int(transformed_data(), meta(), outliers = outliersBX())
      }else{
        boxplot_int_single(transformed_data(), meta(), outliers = outliersBX(), id=id6()) 
      }
    } else if (input$level6 == "Phosphosite"){
      req(data3(), meta2())
      if (meanBX()){
        boxplot_int(data3(), meta2(), outliers = outliersBX())
      }else{
        boxplot_int_single(data3(), meta2(), outliers = outliersBX(), id=id6()) 
      }
    }
  })
  
  output$downloadBoxIntPlot <- downloadHandler(
    filename = function() {
      paste("boxint_plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      width <- input$plotWidth6
      height <- input$plotHeight6
      dpi <- input$plotDPI6
      
      if (meanBX()){
        ggsave(file, plot = boxplot_int(transformed_data(), meta(), outliers = outliersBX()), device = "png", 
               width = width, height = height, units = "cm", dpi = dpi)
      }else{
        ggsave(file, plot = boxplot_int_single(transformed_data(), meta(), outliers = outliersBX()), device = "png", 
               width = width, height = height, units = "cm", dpi = dpi)
      }
    }
  )
  
  observeEvent(input$addText6up, {
    df <- const_df()
    if (input$level6 == "Protein"){
      df <- df[df$Var != "text6upProtein", ]
      new_row <- data.frame(Var = "text6upProtein", Select = input$text6, stringsAsFactors = FALSE)
    } else if (input$level6 == "Phosphosite"){
      df <- df[df$Var != "text6upPhosphosite", ]
      new_row <- data.frame(Var = "text6upPhosphosite", Select = input$text6, stringsAsFactors = FALSE)
    }
    updateTextAreaInput(session, "text6", value = "")
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$addText6down, {
    df <- const_df()
    if (input$level6 == "Protein"){
      df <- df[df$Var != "text6downProtein", ]
      new_row <- data.frame(Var = "text6downProtein", Select = input$text6, stringsAsFactors = FALSE)
    } else if (input$level6 == "Phosphosite"){
      df <- df[df$Var != "text6downPhosphosite", ]
      new_row <- data.frame(Var = "text6downPhosphosite", Select = input$text6, stringsAsFactors = FALSE)
    }
    updateTextAreaInput(session, "text6", value = "")
    const_df(rbind(df, new_row))
  })
  
#### Tab07 - COV Plot ####
  outliersCOV <- reactiveVal(FALSE)
  
  observeEvent(input$toggle_outliersCOV, {
    outliersCOV(!outliersCOV())
  })
  
  output$CovPlot <- renderPlot({
    if (input$level7 == "Protein") {
      req(data(), meta())
      cov_plot(not_transformed_data(), meta(), outliers = outliersCOV())
    } else if (input$level7 == "Phosphosite"){
      req(data3(), meta2())
      cov_plot(not_transformed_data3(), meta2(), outliers = outliersCOV())
    }
  })
  
  output$downloadCovPlot <- downloadHandler(
    filename = function() {
      paste("cov_plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      width <- input$plotWidth7
      height <- input$plotHeight7
      dpi <- input$plotDPI7
      
      ggsave(file, plot = cov_plot(data(), meta(), outliers = outliersCOV()), device = "png", 
             width = width, height = height, units = "cm", dpi = dpi)
    }
  )
  
  observeEvent(input$addText7up, {
    df <- const_df()
    if (input$level7 == "Protein"){
      df <- df[df$Var != "text7upProtein", ]
      new_row <- data.frame(Var = "text7upProtein", Select = input$text7, stringsAsFactors = FALSE)
    } else if (input$level7 == "Phosphosite"){
      df <- df[df$Var != "text7upPhosphosite", ]
      new_row <- data.frame(Var = "text7upPhosphosite", Select = input$text7, stringsAsFactors = FALSE)
    }
    updateTextAreaInput(session, "text7", value = "")
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$addText7down, {
    df <- const_df()
    if (input$level7 == "Protein"){
      df <- df[df$Var != "text7downProtein", ]
      new_row <- data.frame(Var = "text7downProtein", Select = input$text7, stringsAsFactors = FALSE)
    } else if (input$level7 == "Phosphosite"){
      df <- df[df$Var != "text7downPhosphosite", ]
      new_row <- data.frame(Var = "text7downPhosphosite", Select = input$text7, stringsAsFactors = FALSE)
    }
    updateTextAreaInput(session, "text7", value = "")
    const_df(rbind(df, new_row))
  })
  
#### Tab08 - PCA ####
  output$PCAPlot <- renderPlot({
    if (input$level8 == "Protein") {
      req(data(), meta())
      dim_func(data(), meta(), method=input$style8)
    } else if (input$level8 == "Phosphosite"){
      req(data3(), meta2())
      dim_func(data3(), meta2(), method=input$style8)
    }
  })
  
  output$downloadPCAPlot <- downloadHandler(
    filename = function() {
      paste("pca_plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      width <- input$plotWidth8
      height <- input$plotHeight8
      dpi <- input$plotDPI8
      
      ggsave(file, plot = pca_plot(data(), meta()), device = "png", 
             width = width, height = height, units = "cm", dpi = dpi)
    }
  )
  
  observeEvent(input$addText8up, {
    df <- const_df()
    if (input$level8 == "Protein"){
      df <- df[df$Var != "text8upProtein", ]
      new_row <- data.frame(Var = "text8upProtein", Select = input$text8, stringsAsFactors = FALSE)
    } else if (input$level8 == "Phosphosite"){
      df <- df[df$Var != "text8upPhosphosite", ]
      new_row <- data.frame(Var = "text8upPhosphosite", Select = input$text8, stringsAsFactors = FALSE)
    }
    updateTextAreaInput(session, "text8", value = "")
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$addText8down, {
    df <- const_df()
    if (input$level8 == "Protein"){
      df <- df[df$Var != "text8downProtein", ]
      new_row <- data.frame(Var = "text8downProtein", Select = input$text8, stringsAsFactors = FALSE)
    } else if (input$level8 == "Phosphosite"){
      df <- df[df$Var != "text8downPhosphosite", ]
      new_row <- data.frame(Var = "text8downPhosphosite", Select = input$text8, stringsAsFactors = FALSE)
    }
    updateTextAreaInput(session, "text8", value = "")
    const_df(rbind(df, new_row))
  })
  
#### Tab09 - Abundance Plot ####
  observe({
    if (input$level9 == "Protein"){
      req(not_transformed_data())
      updateSelectizeInput(session, "protein9", choices = unique(not_transformed_data()$ProteinNames))
    } else if (input$level9 == "Phosphosite"){
      req(not_transformed_data3())
      updateSelectizeInput(session, "protein9", choices = unique(not_transformed_data3()$PTM_Collapse_key))
    }
  })
  
  observe({
    if (input$level9 == "Protein"){
      req(not_transformed_data())
      updateSelectInput(session, "condition9", 
                        choices = c("All Conditions", unique(meta()$condition)))}
    else if (input$level9 == "Phosphosite"){
      req(data3())
      updateSelectInput(session, "condition9", 
                        choices = c("All Conditions", unique(meta2()$condition)))}
  })
  
  output$abundancePlot <- renderPlotly({
    if (input$level9 == "Protein"){
      req(input$condition9)
      if (input$condition9 == "All Conditions") {
        abundance_plot(not_transformed_data(), meta(), workflow = input$level9)
      } else {
        interactive_abundance_plot(not_transformed_data(), meta(), input$condition9, workflow = input$level9, search = input$protein9)
      }} else if (input$level9 == "Phosphosite"){
        req(input$condition9)
        if (input$condition9 == "All Conditions") {
          abundance_plot(not_transformed_data3(), meta2(), workflow = "Phosphosite")
        } else {
          interactive_abundance_plot(not_transformed_data3(), meta2(), input$condition9, workflow = "Phosphosite", search = input$protein9)
      }}
  })
  
  output$downloadAbPlot <- downloadHandler(
    filename = function() {
      paste("abundance_plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      width <- input$plotWidth9
      height <- input$plotHeight9
      dpi <- input$plotDPI9
      
      if (input$condition9 == "All Conditions") {
        ggsave(file, plot = abundance_plot(data(), meta()), device = "png", 
               width = width, height = height, units = "cm", dpi = dpi)
      } else {
        print("Please work on this!")
      }
    }
  )
  
  observeEvent(input$addText6up, {
    df <- const_df()
    if (input$level6 == "Protein"){
      df <- df[df$Var != "text6upProtein", ]
      new_row <- data.frame(Var = "text6upProtein", Select = input$text6, stringsAsFactors = FALSE)
    } else if (input$level6 == "Phosphosite"){
      df <- df[df$Var != "text6upPhosphosite", ]
      new_row <- data.frame(Var = "text6upPhosphosite", Select = input$text6, stringsAsFactors = FALSE)
    }
    updateTextAreaInput(session, "text6", value = "")
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$addText9down, {
    df <- const_df()
    if (input$level9 == "Protein"){
      df <- df[df$Var != "text9downProtein", ]
      new_row <- data.frame(Var = "text9downProtein", Select = input$text9, stringsAsFactors = FALSE)
    } else if (input$level9 == "Phosphosite"){
      df <- df[df$Var != "text9downPhosphosite", ]
      new_row <- data.frame(Var = "text9downPhosphosite", Select = input$text9, stringsAsFactors = FALSE)
    }
    updateTextAreaInput(session, "text9", value = "")
    const_df(rbind(df, new_row))
  })
  
#### Tab10 - Volcano Plot ####
  observe({
    if (input$level10 == "Protein"){
      req(transformed_data())
      updateSelectInput(session, "condition1_10", 
                        choices = c(unique(meta()$condition)))
    } else if (input$level10=="Phosphosite"){
      req(transformed_data3())
      updateSelectInput(session, "condition1_10", 
                        choices = c(unique(meta2()$condition)))
    }
  })
  
  observe({
    if (input$level10 == "Protein"){
      req(transformed_data())
      updateSelectInput(session, "condition2_10", 
                        choices = c(unique(meta()$condition)))
    } else if (input$level10=="Phosphosite"){
      req(transformed_data3())
      updateSelectInput(session, "condition2_10", 
                        choices = c(unique(meta2()$condition)))
    }
  })
  
  output$VolcPlot <- renderPlotly({
    if (input$level10 == "Protein"){
      req(volcano_data(), meta(), input$condition1_10, input$condition2_10)
      volcano_plot(volcano_data(), meta(), input$condition1_10, input$condition2_10, in_pval = input$in_pval10, in_log2fc = input$in_log2fc10, paired=input$paired10, uncorrected = input$uncorrected10)
    } else if (input$level10 == "Phosphosite"){
      output$VolcPlot <- renderPlotly({
        req(volcano_data3(), meta2(), input$condition1_10, input$condition2_10)
        volcano_plot(volcano_data3(), meta2(), input$condition1_10, input$condition2_10, in_pval = input$in_pval10, in_log2fc = input$in_log2fc10, workflow="Phosphosite", paired=input$paired10, uncorrected = input$uncorrected10)})
    }
  })
  
  data10 <- reactive({
    if (input$level10 == "Protein") {
      req(volcano_data(), meta(), input$condition1_10, input$condition2_10)
      data <- volcano_data_f(volcano_data(), meta(), input$condition1_10, input$condition2_10, in_pval = input$in_pval10, in_log2fc = input$in_log2fc10, paired=input$paired10, uncorrected = input$uncorrected10)
    } else if (input$level10 == "Phosphosite") {
      req(volcano_data3(), meta2(), input$condition1_10, input$condition2_10)
      data <- volcano_data_f(volcano_data3(), meta2(), input$condition1_10, input$condition2_10, in_pval = input$in_pval10, in_log2fc = input$in_log2fc10, workflow="Phosphosite", paired=input$paired10, uncorrected = input$uncorrected10)
    }
    rownames(data) <- NULL
    return(data)
  })
  
  observeEvent(input$addVolc, {
    if(input$level10 == "Protein") {
      req(volcano_data(), meta(), input$condition1_10, input$condition2_10)
      df <- const_df()
      volc_in <- c(input$condition1_10, input$condition2_10, input$in_pval10, input$in_log2fc10)
      volc_in <- paste(volc_in, collapse = "ඞ")
      new_row <- data.frame(Var = "VolcPlot", Select = volc_in, stringsAsFactors = FALSE)
      const_df(rbind(df, new_row))
    } else if (input$level10 == "Phosphosite"){
      req(volcano_data3(), meta2(), input$condition1_10, input$condition2_10)
      df <- const_df()
      phosvolc_in <- c(input$condition1_10, input$condition2_10, input$in_pval10, input$in_log2fc10)
      phosvolc_in <- paste(phosvolc_in, collapse = "ඞ")
      new_row <- data.frame(Var = "PhosVolcPlot", Select = phosvolc_in, stringsAsFactors = FALSE)
      const_df(rbind(df, new_row))
    }
  })
  
  observeEvent(input$addText10up, {
    df <- const_df()
    if (input$level10 == "Protein"){
      df <- df[df$Var != "text10upProtein", ]
      new_row <- data.frame(Var = "text10upProtein", Select = input$text10, stringsAsFactors = FALSE)
    } else if (input$level10 == "Phosphosite"){
      df <- df[df$Var != "text10upPhosphosite", ]
      new_row <- data.frame(Var = "text10upPhosphosite", Select = input$text10, stringsAsFactors = FALSE)
    }
    updateTextAreaInput(session, "text10", value = "")
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$addText10down, {
    df <- const_df()
    if (input$level10 == "Protein"){
      df <- df[df$Var != "text10downProtein", ]
      new_row <- data.frame(Var = "text10downProtein", Select = input$text10, stringsAsFactors = FALSE)
    } else if (input$level10 == "Phosphosite"){
      df <- df[df$Var != "text10downPhosphosite", ]
      new_row <- data.frame(Var = "text10downPhosphosite", Select = input$text10, stringsAsFactors = FALSE)
    }
    updateTextAreaInput(session, "text10", value = "")
    const_df(rbind(df, new_row))
  })
  
  output$table10 <- DT::renderDataTable({
    req(data10())
    DT::datatable(data10(), options = list(pageLength = 20))
  })
  
  output$download10<- downloadHandler(
    filename = function() {
      paste("volcanodata_", Sys.Date(), "_", input$condition1_10 ,"_", input$condition2_10 ,".csv", sep = "")
    },
    content = function(file) {
      write.csv(data10(), file, row.names = FALSE)
    }
  )
  
#### Tab11 - GSEA ####
  different_genes_df <- reactive({
    req(volcano_data(), meta(), input$condition1_10, input$condition2_10)
    different_genes(volcano_data(), meta(), input$condition1_10, input$condition2_10, in_pval = input$in_pval10, in_log2fc = input$in_log2fc10)
  })
  
  output$UpregEnrichmentPlot <- renderPlot({
    req(different_genes_df())
    enrichment_analysis(different_genes_df()$Upregulated, top_n = input$top_n11, min_num = input$filter11min, max_num = input$filter11max)
  })
  
  output$DownregEnrichmentPlot <- renderPlot({
    req(different_genes_df())
    enrichment_analysis(different_genes_df()$Downregulated, top_n = input$top_n11, min_num = input$filter11min, max_num = input$filter11max)
  })
  
  output$downloadUpPlot <- downloadHandler(
    filename = function() {
      paste("enrichup_plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      width <- input$plotWidth11
      height <- input$plotHeight11
      dpi <- input$plotDPI11
      
      ggsave(file, plot = enrichment_analysis(different_genes_df()$Upregulated, top_n = input$top_n11, min_num = input$filter11min, max_num = input$filter11max), device = "png", 
             width = width, height = height, units = "cm", dpi = dpi)
    }
  )
  
  output$downloadDownPlot <- downloadHandler(
    filename = function() {
      paste("enrichdown_plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      width <- input$plotWidth11
      height <- input$plotHeight11
      dpi <- input$plotDPI11
      
      ggsave(file, plot = enrichment_analysis(different_genes_df()$Downregulated, top_n = input$top_n11, min_num = input$filter11min, max_num = input$filter11max), device = "png", 
             width = width, height = height, units = "cm", dpi = dpi)
    }
  )
  
  output$table11up <- DT::renderDataTable({
    req(different_genes_df())
    DT::datatable(different_genes_df()$Upregulated, options = list(pageLength = 5))
  })
  
  output$table11down <- DT::renderDataTable({
    req(different_genes_df())
    DT::datatable(different_genes_df()$Downregulated, options = list(pageLength = 5))
  })
  
#### Tab12 - Correlation Plot ####
  id12 <- reactiveVal(TRUE)
  
  observeEvent(input$toggle_id12, {
    id12(!id12())
  })
  
  MeCorr <- reactiveVal(FALSE)
  
  observeEvent(input$Change12, {
    MeCorr(!MeCorr())
  })
  
  output$CorrPlot <- renderPlot({
    if (input$level12 == "Protein") {
      req(data(), meta())
      corr_plot(data(), meta(), MeCorr(), id12())
    } else if (input$level12 == "Phosphosite"){
      req(data3(), meta2())
      corr_plot(data3(), meta2(), MeCorr(), id12())
    }
  })
  
  output$downloadCorrPlot <- downloadHandler(
    filename = function() {
      paste("corr_plot", Sys.Date(), ".png", sep = "")
    },
    
    content = function(file) {
      width <- input$plotWidth12 / 2.54 
      height <- input$plotHeight12 / 2.54 
      dpi <- input$plotDPI12
      
      png(file, width = width, height = height, units = "in", res = dpi)
      corr_plot(data(), meta(), MeCorr(), id12())
      dev.off()
    }
  )
  
  observeEvent(input$addText12up, {
    df <- const_df()
    if (input$level12 == "Protein"){
      df <- df[df$Var != "text12upProtein", ]
      new_row <- data.frame(Var = "text12upProtein", Select = input$text12, stringsAsFactors = FALSE)
    } else if (input$level12 == "Phosphosite"){
      df <- df[df$Var != "text12upPhosphosite", ]
      new_row <- data.frame(Var = "text12upPhosphosite", Select = input$text12, stringsAsFactors = FALSE)
    }
    updateTextAreaInput(session, "text12", value = "")
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$addText12down, {
    df <- const_df()
    if (input$level12 == "Protein"){
      df <- df[df$Var != "text12downProtein", ]
      new_row <- data.frame(Var = "text12downProtein", Select = input$text12, stringsAsFactors = FALSE)
    } else if (input$level12 == "Phosphosite"){
      df <- df[df$Var != "text12downPhosphosite", ]
      new_row <- data.frame(Var = "text12downPhosphosite", Select = input$text12, stringsAsFactors = FALSE)
    }
    updateTextAreaInput(session, "text12", value = "")
    const_df(rbind(df, new_row))
  })
  
#### Tab13 - Heatmap ####
  id13 <- reactiveVal(TRUE)
  show13 <- reactiveVal(FALSE)
  
  observeEvent(input$toggle_id13, {
    id13(!id13())
  })
  
  observeEvent(input$show13, {
    show13(!show13())
  })
  
  observeEvent(input$fix13, {
    clear_all_plots()
  })
  
  output$heatPlot <- renderPlot({
    if (input$level13=="Protein"){
      req(data(), meta())
      if (show13()==TRUE){
        heatmap_plot(data(), meta(), id13())
      } else if (show13()==FALSE){
        heatmap_plot_nmv(transformed_data(), meta(), id13())
        }
    } else if (input$level13=="Phosphosite"){
      req(data3(), meta2())
      if (show13()==TRUE){
        heatmap_plot(data3(), meta2(), id13())
      } else if (show13()==FALSE){
        heatmap_plot_nmv(transformed_data3(), meta2(), id13())
      }
    }
  })
  
  output$downloadheatPlot <- downloadHandler(
    filename = function() {
      paste("heatPlot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      width <- input$plotWidth13
      height <- input$plotHeight13 / 2.54 
      dpi <- input$plotDPI13 / 2.54 
      
      png(file, width = width, height = height, units = "in", res = dpi)
      heatmap_plot(data(), meta(), id13())
      dev.off()
    }
  )
  
  observeEvent(input$addText13up, {
    df <- const_df()
    if (input$level13 == "Protein"){
      df <- df[df$Var != "text13upProtein", ]
      new_row <- data.frame(Var = "text13upProtein", Select = input$text13, stringsAsFactors = FALSE)
    } else if (input$level13 == "Phosphosite"){
      df <- df[df$Var != "text13upPhosphosite", ]
      new_row <- data.frame(Var = "text13upPhosphosite", Select = input$text13, stringsAsFactors = FALSE)
    }
    updateTextAreaInput(session, "text13", value = "")
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$addText13down, {
    df <- const_df()
    if (input$level13 == "Protein"){
      df <- df[df$Var != "text13downProtein", ]
      new_row <- data.frame(Var = "text13downProtein", Select = input$text13, stringsAsFactors = FALSE)
    } else if (input$level13 == "Phosphosite"){
      df <- df[df$Var != "text13downPhosphosite", ]
      new_row <- data.frame(Var = "text13downPhosphosite", Select = input$text13, stringsAsFactors = FALSE)
    }
    updateTextAreaInput(session, "text13", value = "")
    const_df(rbind(df, new_row))
  })
  
#### Tab14 - RT Plot ####
  line14 <- reactiveVal(F)
  
  observeEvent(input$line14, {
    line14(!line14())
  })
  
  output$RTPlot <- renderPlot({
    req(data2(), meta())
    RTvspredRT_plot(data2(), meta(), method=input$style14, add_line=line14(), bin = input$bins14)
  })
  
  observeEvent(input$addText14up, {
    df <- const_df()
    df <- df[df$Var != "text14up", ]
    new_row <- data.frame(Var = "text14up", Select = input$text14, stringsAsFactors = FALSE)
    updateTextAreaInput(session, "text14", value = "")
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$addText14down, {
    df <- const_df()
    df <- df[df$Var != "text14down", ]
    new_row <- data.frame(Var = "text14down", Select = input$text14, stringsAsFactors = FALSE)
    updateTextAreaInput(session, "text14", value = "")
    const_df(rbind(df, new_row))
  })
  
#### Tab15 - Modification Plot ####
  id15 <- reactiveVal(TRUE)
  
  observeEvent(input$toggle_id15, {
    id15(!id15())
  })
  
  output$ModPlot <- renderPlot({
    req(data2(), meta())
    modification_plot(data2(), meta(), id15())
  })
  
  observeEvent(input$addText15up, {
    df <- const_df()
    df <- df[df$Var != "text15up", ]
    new_row <- data.frame(Var = "text15up", Select = input$text15, stringsAsFactors = FALSE)
    updateTextAreaInput(session, "text15", value = "")
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$addText15down, {
    df <- const_df()
    df <- df[df$Var != "text15down", ]
    new_row <- data.frame(Var = "text15down", Select = input$text15, stringsAsFactors = FALSE)
    updateTextAreaInput(session, "text15", value = "")
    const_df(rbind(df, new_row))
  })
  
#### Tab16 - Missed Cleavage Plot ####
  id16 <- reactiveVal(TRUE)
  
  observeEvent(input$toggle_id16, {
    id16(!id16())
  })
  
  output$MCPlot <- renderPlot({
    req(data2(), meta())
    missed_cl_plot(data2(), meta(), id16())
  })
  
  observeEvent(input$addText16up, {
    df <- const_df()
    df <- df[df$Var != "text16up", ]
    new_row <- data.frame(Var = "text16up", Select = input$text16, stringsAsFactors = FALSE)
    updateTextAreaInput(session, "text16", value = "")
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$addText16down, {
    df <- const_df()
    df <- df[df$Var != "text16down", ]
    new_row <- data.frame(Var = "text16down", Select = input$text16, stringsAsFactors = FALSE)
    updateTextAreaInput(session, "text16", value = "")
    const_df(rbind(df, new_row))
  })
  
#### Tab17 - Boxplot, single protein ####
  observe({
    if (input$level17 == "Protein"){
      req(transformed_data())
      updateSelectizeInput(session, "protein17", choices = unique(data()$ProteinNames))
      updateSelectizeInput(session, "conditions17", choices = unique(meta()$condition))
    } else if (input$level17 == "Phosphosite"){
      req(transformed_data())
      updateSelectizeInput(session, "protein17", choices = unique(data3()$PTM_Collapse_key))
      updateSelectizeInput(session, "conditions17", choices = unique(meta2()$condition))
    }
  })
  
  output$ProtBoxPlot <- renderPlot({
    if (input$level17 == "Protein"){
      req(transformed_data(), meta())
      validate(
        need(length(input$conditions17) >= 2, "Please select at least two conditions")
      )
      plot <- tryCatch({
        compare_prot_box(transformed_data(), meta(), input$conditions17, input$protein17)
      }, error = function(e) {
        stop(safeError(e))
      })
      print(plot)
    } else if (input$level17=="Phosphosite"){
      req(transformed_data3(), meta2())
      validate(
        need(length(input$conditions17) >= 2, "Please select at least two conditions")
      )
      plot <- tryCatch({
        compare_prot_box(transformed_data3(), meta2(), input$conditions17, input$protein17, workflow = "Phosphosite")
      }, error = function(e) {
        stop(safeError(e))
      })
      print(plot)
    }
  })
  
#### Tab18 - Lineplot, single protein ####
  observe({
    if (input$level18 == "Protein"){
      req(transformed_data())
      updateSelectizeInput(session, "protein18", choices = unique(data()$ProteinNames))
      updateSelectizeInput(session, "conditions18", choices = unique(meta()$condition))
    } else if (input$level18 == "Phosphosite"){
      req(transformed_data())
      updateSelectizeInput(session, "protein18", choices = unique(data3()$PTM_Collapse_key))
      updateSelectizeInput(session, "conditions18", choices = unique(meta2()$condition))
    }
  })
  
  output$ProtLinePlot <- renderPlot({
    if (input$level18=="Protein"){
      req(transformed_data(), meta())
      plot <- tryCatch({
        compare_prot_line(transformed_data(), meta(), input$conditions18, input$protein18)
      }, error = function(e) {
        stop(safeError(e))
      })
      print(plot)
    } else if (input$level18=="Phosphosite"){
      req(transformed_data3(), meta2())
      plot <- tryCatch({
        compare_prot_line(transformed_data3(), meta2(), input$conditions18, input$protein18, workflow = "Phosphosite")
      }, error = function(e) {
        stop(safeError(e))
      })
      print(plot)
    }
  })
  
#### Tab19 - Phossite Plot ####
  output$PhossitePlot <- renderPlot({
    req(data3())
    simple_phos_site_plot(data3(), filter=input$cutoff19)
  })
  
  observeEvent(input$addText19up, {
    df <- const_df()
    df <- df[df$Var != "text19up", ]
    new_row <- data.frame(Var = "text19up", Select = input$text19, stringsAsFactors = FALSE)
    updateTextAreaInput(session, "text19", value = "")
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$addText19down, {
    df <- const_df()
    df <- df[df$Var != "text19down", ]
    new_row <- data.frame(Var = "text19down", Select = input$text19, stringsAsFactors = FALSE)
    updateTextAreaInput(session, "text19", value = "")
    const_df(rbind(df, new_row))
  })
  
#### Tab20.1 - KSEA ####
  observe({
    req(transformed_data3())
    updateSelectInput(session, "condition1_20", 
                      choices = c(unique(meta2()$condition)))
  })
  
  observe({
    req(transformed_data3())
    updateSelectInput(session, "condition2_20", 
                      choices = c(unique(meta2()$condition)))
  })
  
  output$downKSEA <- downloadHandler(
    filename = function() {
      paste("KSEA_data_", input$condition1_20, "_", input$condition2_20, ".txt", sep = "")  
    },
    content = function(file) {
      write.table(prepare_KSEA(transformed_data3(), meta2(), input$condition1_20, input$condition2_20), file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )
  
  KSEA_data <- reactive({
    req(input$KSEA_data)
    file <- input$KSEA_data$datapath
    read.delim(file, sep = "\t")
  })
  
  output$KinaseVolcPlot <- renderPlotly({
    req(KSEA_data(), input$in_pval20)
    kinase_volcano(KSEA_data(), input$in_pval20)
  })
  
  output$downTREE <- downloadHandler(
    filename = function() {
      paste("TREE_data_", input$condition1_20, "_", input$condition2_20, ".zip", sep = "")  
    },
    content = function(file) {
      temp_dir <- tempdir()
      
      file1 <- file.path(temp_dir, "enrichment_scores.txt")
      file2 <- file.path(temp_dir, "pvals.txt")
      
      write.table(prepare_tree_data_es(KSEA_data()), 
                  file1, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
      
      write.table(prepare_tree_data_pv(KSEA_data()), 
                  file2, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
      
      zip::zipr(file, files = c(file1, file2))
    }
  )
  
  output$treePlot <- renderUI({
    req(input$TreeSVG)
    svg_content <- readLines(input$TreeSVG$datapath)
    HTML(paste(svg_content, collapse = "\n"))
  })
  
#### Tab20.2 - KSEA (Kinact)####
  observe({
    req(transformed_data3())
    updateSelectInput(session, "condition1_202", 
                      choices = c(unique(meta2()$condition)))
  })
  
  observe({
    req(transformed_data3())
    updateSelectInput(session, "condition2_202", 
                      choices = c(unique(meta2()$condition)))
  })
  
  data202 <- reactive({
    req(volcano_data3(), meta2(), input$condition1_202, input$condition2_202)
    volc_data <- volcano_data_f(volcano_data3(), meta2(), input$condition1_202, input$condition2_202, workflow="Phosphosite")
    rownames(volc_data) <- NULL
    org_data = volcano_data3()
    data = match_volc_and_org_data(volc_data = volc_data, org_data = org_data)
    return(data)
  })
  
  output$KinactPlot <- renderPlotly({
      req(data202())
      kinact_kinase_activity(data202(), top_n = input$top_n202, NetworKIN = input$NetworKIN202, NetworKIN.cutoff = input$NetworKIN.cutoff202, m.cutoff = input$m.cutoff202)
    })
  
  output$KinactPlotSites <- renderPlotly({
    req(data202(), input$Kinase_202)
    downstream_phossite_volc(data202(), input$Kinase_202, NetworKIN = input$NetworKIN202)
  })
  
#### Tab21 - Phossite Coverage Plot ####
  output$PhossiteCoveragePlot <- renderPlot({
    req(data3(), meta2())
    phossite_coverage_plot(data3(), meta2())
  })
  
#### Tab22 - Protein centric ####
  db_fasta <- reactiveVal(NULL)
  
  observeEvent(input$db_fasta22, {
    #db_fasta_stat <- readFASTA("C:/Users/jonas/OneDrive/Desktop/Database/Uniprot/human_prot_uniprotkb_proteome_UP000005640_2024_04_11.fasta") #Human
    db_fasta_stat <- readFASTA("C:/Users/jonas/OneDrive/Desktop/Database/Uniprot/UP000000589_10090.fasta") #Mouse
    db_fasta_stat = as.data.frame(db_fasta_stat)
    db_fasta_stat = as.data.frame(t(db_fasta_stat))
    db_fasta_stat$name <- rownames(db_fasta_stat)
    db_fasta_stat$name <- sapply(strsplit(as.character(db_fasta_stat$name), "\\."), tail, 1)
    db_fasta(db_fasta_stat)
  })
  
  observe({
    req(transformed_data())
    updateSelectizeInput(session, "protein20", choices = unique(data()$ProteinNames))
  })
  
  output$seq_allign <- renderPrint({
    split_string <- strsplit(input$protein20, ";")[[1]]
    vis_coverage(data2(), split_string[1], chunk_size = input$chunk20, db = db_fasta())
  })
  
  output$plot_3d22 <- renderR3dmol({
    protein <- input$protein20
    plot <- model_3d(data2(), protein, db = db_fasta())
    return(plot)
  })
  
#### Tabn1 - Impute data ####
  output$Impute1 <- renderPlot({
    if (input$leveln1 == "Protein") {
      req(transformed_data(), meta())
      impute_values(transformed_data(), meta(), ret=1, q=input$qn1, adj_std = input$adj_stdn1, seed=input$seedn1)
    } else if (input$leveln1 == "Phosphosite"){
      req(transformed_data3(), meta2())
      impute_values(transformed_data3(), meta2(), ret=1, q=input$qn1, adj_std = input$adj_stdn1, seed=input$seedn1)
    }
  })
  
  output$Impute2 <- renderPlot({
    if (input$leveln1 == "Protein") {
      req(transformed_data(), meta())
      impute_values(transformed_data(), meta(), ret=2, q=input$qn1, adj_std = input$adj_stdn1, seed=input$seedn1)
    } else if (input$leveln1 == "Phosphosite"){
      req(transformed_data3(), meta2())
      impute_values(transformed_data3(), meta2(), ret=2, q=input$qn1, adj_std = input$adj_stdn1, seed=input$seedn1)
    }
  })
  
  output$Impute3 <- renderPlot({
    if (input$leveln1 == "Protein") {
      req(transformed_data(), meta())
      impute_values(transformed_data(), meta(), ret=3, q=input$qn1, adj_std = input$adj_stdn1, seed=input$seedn1)
    } else if (input$leveln1 == "Phosphosite"){
      req(transformed_data3(), meta2())
      impute_values(transformed_data3(), meta2(), ret=3, q=input$qn1, adj_std = input$adj_stdn1, seed=input$seedn1)
    }
  })
  
  imputed_data <- reactiveVal(NULL)
  imputed_data3 <- reactiveVal(NULL)
  
  observeEvent(input$ImputeEVE, {
    if (input$leveln1 == "Protein") {
      req(transformed_data(), meta())
      imputed_data(impute_values(transformed_data(), meta(), ret=0, q=input$qn1, adj_std = input$adj_stdn1, seed=input$seedn1))
      volcano_data <- reactive({
        imputed_data()
      })
    } else if (input$leveln1 == "Phosphosite"){
      req(transformed_data3(), meta2())
      imputed_data3(impute_values(transformed_data3(), meta2(), ret=0, q=input$qn1, adj_std = input$adj_stdn1, seed=input$seedn1))
      volcano_data3 <- reactive({
        imputed_data3()
      })
    }
  })
  
  output$imputed_data_tab <- DT::renderDataTable({
    if (input$leveln1 == "Protein") {
      req(imputed_data())
      DT::datatable(imputed_data(), options = list(pageLength = 20))
    } else if (input$leveln1 == "Phosphosite"){
      req(imputed_data3())
      DT::datatable(imputed_data3(), options = list(pageLength = 20))
    }
  })
  
  output$imputed_data_down <- downloadHandler(
      filename = function() {
        if (input$leveln1 == "Protein") {
          paste("imputed_data-", Sys.Date(), ".csv", sep = "")
        } else if (input$leveln1 == "Phosphosite") {
          paste("imputed_data_phospho-", Sys.Date(), ".csv", sep = "")
        }
      },
      content = function(file) {
        if (input$leveln1 == "Protein") {
          req(imputed_data())
          write.csv(imputed_data(), file, row.names = FALSE)
        } else if (input$leveln1 == "Phosphosite") {
          req(imputed_data3())
          write.csv(imputed_data3(), file, row.names = FALSE)
        }
      }
  )
  
#### TabSummary - Summary ####
  output$transformed_data_prot <- DT::renderDataTable({
    req(transformed_data())
    DT::datatable(transformed_data(), options = list(pageLength = 5))
  })
  
  output$meta_object <- DT::renderDataTable({
    req(meta())
    DT::datatable(meta(), options = list(pageLength = 5))
  })
  
  output$collapse_data_sum <- DT::renderDataTable({
    req(data3())
    DT::datatable(data3(), options = list(pageLength = 5))
  })
  
  output$meta_object2 <- DT::renderDataTable({
    req(meta2())
    DT::datatable(meta2(), options = list(pageLength = 5))
  })
  
  output$download_meta <- downloadHandler(
    filename = function() {
      paste("metadata-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(meta(), file, row.names = FALSE)
    }
  )
  
  output$download_meta2 <- downloadHandler(
    filename = function() {
      paste("metadata_phos-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(meta2(), file, row.names = FALSE)
    }
  )
  
  output$download_collapse_data_sum <- downloadHandler(
    filename = function() {
      paste("collapse_data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(data3(), file, row.names = FALSE)
    }
  )
  
#### LogLogic ####
  rv <- reactiveValues(data = NULL)
  
  observe({
    var_list <- list("Transformed", "Transformed3",
                     "CoveragePlotID", "BoxplotIntID", "BoxplotIntOut", "BoxplotIntMean",
                     "CovPlotOut", "CorrPlotDisplay", "CorrPlotID", "HeatmapID",
                     "RTPlotSytle", "ModPlotID", "MissedCleavID", "HexbinsRT", "RTline")
    
    select_list <- list(was_transformed(), was_transformed3(),
                        id3(), id6(), outliersBX(), 
                        meanBX(), outliersCOV(), MeCorr(), 
                        id12(), id13(), input$style14, 
                        id15(), id16(), input$bins14, line14())
    
    df <- do.call(rbind, Map(data.frame, Var = var_list, Select = select_list, stringsAsFactors = FALSE))
    df$Select[df$Select == 1] <- TRUE
    df$Select[df$Select == 0] <- FALSE
    rv$data <- df
  })
  
  const_df <- reactiveVal(data.frame(Var = character(), Select = character(), stringsAsFactors = FALSE))
  
  observeEvent(input$file, {
    req(input$file)
    df <- const_df()
    file_name <- input$file$name
    df <- df[df$Var != "Data1", ]
    new_row <- data.frame(Var = "Data1", Select = file_name, stringsAsFactors = FALSE)
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$file2, {
    req(input$file2)
    df <- const_df()
    file_name <- input$file2$name
    df <- df[df$Var != "Data2", ]
    new_row <- data.frame(Var = "Data2", Select = file_name, stringsAsFactors = FALSE)
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$file3, {
    req(input$file3)
    df <- const_df()
    file_name <- input$file3$name
    df <- df[df$Var != "Data3", ]
    new_row <- data.frame(Var = "Data3", Select = file_name, stringsAsFactors = FALSE)
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$KSEA_data, {
    req(input$KSEA_data)
    df <- const_df()
    file_name <- input$KSEA_data$name
    new_row <- data.frame(Var = "KSEA_data", Select = file_name, stringsAsFactors = FALSE)
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$TreeSVG, {
    req(input$TreeSVG)
    df <- const_df()
    file_name <- input$TreeSVG$name
    new_row <- data.frame(Var = "tree_path", Select = file_name, stringsAsFactors = FALSE)
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$upload_meta, {
    req(input$upload_meta)
    df <- const_df()
    file_name <- input$upload_meta$name
    df <- df[df$Var != "Meta1", ]
    new_row <- data.frame(Var = "Meta1", Select = file_name, stringsAsFactors = FALSE)
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$upload_meta2, {
    req(input$upload_meta2)
    df <- const_df()
    file_name <- input$upload_meta2$name
    df <- df[df$Var != "Meta2", ]
    new_row <- data.frame(Var = "Meta2", Select = file_name, stringsAsFactors = FALSE)
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$apply_filter, {
    req(data(), meta(), input$filter_num)
    df <- const_df()
    filter_num <- input$filter_num
    df <- df[df$Var != "FilterNum", ]
    new_row <- data.frame(Var = "FilterNum", Select = filter_num, stringsAsFactors = FALSE)
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$apply_filter2, {
    req(data3(), meta2(), input$filter_num2)
    df <- const_df()
    filter_num2 <- input$filter_num2
    df <- df[df$Var != "FilterNum2", ]
    new_row <- data.frame(Var = "FilterNum2", Select = filter_num2, stringsAsFactors = FALSE)
    const_df(rbind(df, new_row))
  })
  
  combined_df <- reactive({
    req(rv$data)
    if (nrow(const_df()) == 0) {
      return(rv$data)
    }
    rbind(rv$data, const_df())
  })
  
  output$log <- downloadHandler(
    filename = function() {
      paste("report-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(combined_df(), file, row.names = FALSE)
    }
  )
  
  output$LogTable <- renderTable({
    df <- combined_df() 
    if (is.data.frame(df)) {  
      xtable(df)  
    } else {
      NULL 
    }
  })
}

#### ShinyApp logic ####
shinyApp(ui = ui, server = server)