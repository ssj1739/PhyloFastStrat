# Phylostratigraphy App
VERSION = "0.3.2a"

# Requisite libraries:
library(BiocManager)
options(repos = BiocManager::repositories())
library(shiny)
library(shinyjs)
library(shinycssloaders)
library(shinyWidgets)
library(queryup)
library(readxl)
library(org.Hs.eg.db)
library(phylostratr)
library(tidyverse)
library(topGO)
library(GO.db)
library(fgsea)
library(AnnotationDbi)
library(pbapply)
library(STRINGdb)
library(quarto)
library(DT)
library(tidyverse)
library(shinythemes)
library(googlesheets4)
library(ggtext)


##### PREPROCESSING: #####
# Load pre-processed data:
# Gene list:
# MasterGeneLists <- readRDS("data/MasterGeneLists_Original.rds")
# if(file.exists("data/MasterGeneLists_Updated.rds")){
#   MasterGeneLists <- readRDS("data/MasterGeneLists_Updated.rds")
# }
googlesheets4::gs4_auth(cache = F)
MasterGeneLists <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1oEktKQNY8kaEFpImEW-wnoBz4dqLa3eWgrji23AFU4U/edit?gid=1728211307#gid=1728211307")


# strata <- readRDS("9606_strata.rds")
# Merge results into a single hittable
#results <- merge_besthits(strata)
#ph <- stratify(results, classify_by_adjusted_pvalue(0.001))
#saveRDS(ph, file = "ph_0.001.rds")
ph <- readRDS("data/ph_0.001.rds")

levels(ph$mrca_name)[match("cellular organisms", levels(ph$mrca_name))] <- "Cellular Organisms"

Final_Taxa_Emergence_Timeline_with_Cellular_Organisms <- read.csv("data/Final_Taxa_Emergence_Timeline_with_Cellular_Organisms.csv")

levels(ph$mrca_name) <- paste0(levels(ph$mrca_name), ": ", paste0(Final_Taxa_Emergence_Timeline_with_Cellular_Organisms$Emergence..MYA., " MYA"))

#pathways <- fgsea::gmtPathways("c5.go.v2023.2.Hs.symbols.gmt")
#saveRDS(pathways, file = "PhylostratigraphyApp/c5.go.v2023.2.Hs.symbols.rds")
pathways <- readRDS("data/c5.go.v2023.2.Hs.symbols.rds")

# up <- UniProt.ws::UniProt.ws()
# gene_uniprot_conv_df <- UniProt.ws::select(up, keys = ph$qseqid, to = "Gene_Name", 
#                                            columns = c("Entry", "Gene_Name"))
# colnames(gene_uniprot_conv_df) <- c("Entry", "Gene")

#saveRDS(gene_uniprot_conv_df, "data/uniprot_gene_conv_df_updated.rds")
#gene_uniprot_conv_df <- readRDS("PhylostratigraphyApp/uniprot_gene_conv_df.rds")
#uniprotkb <- read_delim("../uniprotkb_organism_id_9606_2025_03_15.tsv.gz", 
#                                                                    delim = "\t", escape_double = FALSE,
#                                                                    trim_ws = TRUE)

#gene_uniprot_conv_df <- merge(uniprotkb, gene_uniprot_conv_df, by = "Entry", all = T)
#gene_uniprot_conv_df$Gene <- sapply(gene_uniprot_conv_df$`Gene Names`, function(x) strsplit(x, split = " ")[[1]][1])
#saveRDS(gene_uniprot_conv_df, file = "data/uniprot_gene_conv_df_updated.rds")
gene_uniprot_conv_df <- readRDS("data/uniprot_gene_conv_df_updated.rds")
#genes <- AnnotationDbi::select(org.Hs.eg.db, gene_uniprot_conv_df$`Entry`, "SYMBOL","UNIPROT")
#gene_uniprot_conv_df$Gene <- genes$SYMBOL[match(gene_uniprot_conv_df$Entry, genes$UNIPROT)]
#saveRDS(gene_uniprot_conv_df, file = "data/uniprot_gene_conv_df_updated.rds")

ph <- merge(ph, gene_uniprot_conv_df, by.x = "qseqid", by.y = "Entry", all.x = T)
ph$Gene <- ph$`Gene Names (primary)`
ph$gene <- ph$Gene

all_genes <- ph$Gene[!is.na(ph$Gene)]
names(all_genes) <- NULL
all_genes <- unique(all_genes)

ph_filt <- ph |>
  dplyr::filter(!duplicated(Gene))

my_theme <- ggplot2::theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 12),
    plot.title = element_text(face = "bold", size = 14)
  )

##### UI #####
ui <- fluidPage(
  theme = shinytheme('yeti'),
  # Application title
  titlePanel(paste0("PhyloFastStrat v", VERSION)),
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      shinyWidgets::switchInput(inputId = "runmode", label = "Run Mode", offLabel = "Preset Gene Lists", onLabel = "Custom Gene Lists", value = FALSE),
      conditionalPanel(
        condition = "!input.runmode",
        selectInput(
          inputId = "disease_of_interest",
          label = "Disease-Gene List:",
          choices = c(colnames(MasterGeneLists)),
          multiple = F,
          selected = NULL,
          selectize = T
        )
      ),
      conditionalPanel(
        condition = "input.runmode",
        textInput(inputId = "other_disease_of_interest", "Gene set name:", value = "My Gene List", placeholder = "Gene List"),
        selectizeInput(inputId = "selectized_genes", label = "Genes to query:",
                       choices = NULL,
                       multiple = TRUE,
                       selected = NULL,
                       options = list(
                         splitOn = I("(function() { return/[,; ]/; })()"),
                         plugins = list('remove_button'),
                         delimiter = " ",
                         persist = TRUE
                       )),
        fileInput(
          inputId = "custom_file",
          label = "Custom Gene List:",
          accept = ".csv",
          buttonLabel = "Upload"
        ),
        # a(href="genelist.csv", "Download example gene list.", download=NA, target="_blank"),
        # tags$br(),
        # tags$br(),
        #actionButton(inputId = "saverds", label = "Cache custom data"),
        #textOutput("cache_message")
      ),
      headerPanel(""),
      htmlOutput("numGenesMappedUP"),
      headerPanel(""),
      shinyWidgets::switchInput(inputId = "no_isoforms", label = "Exclude alternative isoforms", onLabel = "Yes", offLabel = "No", value = FALSE),
      htmlOutput("numIsoformsMappedUP"),
      headerPanel(""),
      actionButton("run", "Run analysis on selected genes"),

    ),
    
    # Show a plot of the generated distribution
    mainPanel(tabsetPanel(
      type = "tabs",
      tabPanel("Origin of genes",
               plotOutput("ph_plot", click = "ph_plot_click"),
               DT::dataTableOutput("selected_ph", height = "100%"),
               ),
      tabPanel("Rate comparison",
               plotOutput("rate_plot"),
               "Rates are calculated by dividing the number of genes in a given phylostrata by the total number of genes in that phylostrata.\n 
    The odds ratio is calculated by dividing the rate of the disease gene set by the rate of the overall human genome.\n 
    A value greater than 1 indicates that the disease gene set is enriched in that phylostrata.\n
    The presence of an asterisk ('*') indicates that the odds ratio is significantly different from 1 (p < 0.05)."
               ),
      tabPanel("Results table",
               shinyWidgets::switchInput(inputId = "show_uniprot_ids", 
                                         label = "Show individual isoform IDs?", 
                                         onLabel = "Yes", 
                                         offLabel = "No", 
                                         value = FALSE),
               DT::dataTableOutput("phylo_table", height = "100%")
               ),
      #tabPanel("Gene map",
       #        selectInput("ps", label = "Select a phylostrata to visualize:", choices = unique(ph$mrca_name)),
               # actionButton("retrieveSTRING", "Retrieve STRING plot", value = NULL),
               # useShinyjs(),
               # extendShinyjs(script = "www/my_functions.js", functions = "loadStringData"),
               # includeScript("http://string-db.org/javascript/combined_embedded_network_v2.0.2.js"),
               # includeScript("https://blueimp.github.io/JavaScript-Canvas-to-Blob/js/canvas-to-blob.js"),
               # includeScript("https://cdnjs.cloudflare.com/ajax/libs/canvg/1.4/rgbcolor.min.js"),
               # includeScript("https://cdnjs.cloudflare.com/ajax/libs/stackblur-canvas/1.4.1/stackblur.min.js"),
               # includeScript("https://cdn.jsdelivr.net/npm/canvg/dist/browser/canvg.min.js"),
               # includeCSS("www/style.css"),
               #plotOutput("string_plot", height = "800px") %>% withSpinner()
       #        ),
      tabPanel(
        "Functional enrichment plots",
        #plotly::plotlyOutput("GOplot2", width = "100%", height = "700px") %>% withSpinner(),
        numericInput(inputId = "labeltop", label = "Number of top GO-BPs to label:", value = 5, min = 1, max = 10),
        numericInput(inputId = "cappval", label = "Max -log P-value to visualize", value = 100, min = 1, max = 100),
        sliderInput(inputId = "pval_range", label = "P-value range to visualize", min = 0, max = 1, value = c(0, 0.05), step = 0.001, pre = "P-value range: ", ),
        plotOutput("GOplot_static", height = "700px", click = "go_plot_click") %>% withSpinner(),
        plotOutput("selected_go", height = "700px") %>% withSpinner(),
        #plotly::plotlyOutput("duplicated_go_plot", width = "100%", height = "500px") %>% withSpinner()
      )
    )
    )
  )
)

##### SERVER #####
server <- function(input, output, session) {
  # Dynamic selectize on server-side
  updateSelectizeInput(session, inputId = "selectized_genes", choices = all_genes, server = TRUE, selected = NULL)

  # If run mode is switched, change disease of interest selection
  # observeEvent(input$runmode, {
  #   updateSelectInput(session, "disease_of_interest")
  #   #print("Triggered reset.")
  #   updateTextInput(session, "mytext", value = "test")
  # })
  
  # Set genes to query based on input, and refresh with "run" button
  ##### genes_of_interest ####
  genes_of_interest_reactive <- eventReactive(input$run, {
    if(!is.null(input$selectized_genes) & input$runmode){
      genes_of_interest <- input$selectized_genes
    }else if(!is.null(input$disease_of_interest) & !input$runmode){
      genes_of_interest <- MasterGeneLists[[input$disease_of_interest]]
    }else{
      custom_gene_list <- readr::read_csv(file = input$custom_file$datapath)
      genes_of_interest <- custom_gene_list[,1]
    }
    # Save a copy of genes_of_interest to compare input to what is processed
    genes_of_interest_input <- genes_of_interest[!is.na(genes_of_interest)]
    
    # Clean genes of interest
    genes_of_interest <-
      genes_of_interest[!is.na(genes_of_interest)]
    genes_of_interest <-
      genes_of_interest[!grepl("RNU", genes_of_interest)]
    genes_of_interest <- unique(genes_of_interest)
    genes_of_interest <-
      sapply(genes_of_interest, function(x)
        strsplit(x, split = " ", fixed = T)[[1]][1])
    
    # How many genes were found in the phylostrata object?
    genes_of_interest_found <- unique(ph$Gene[ph$Gene %in% genes_of_interest])
    
    ##### Main mapping messages #####
    output$numGenesMappedUP <- renderUI({
      # Message to show gene mapping to phylostrat results
      HTML(paste0(length(genes_of_interest_input), " input genes; ", length(genes_of_interest), " gene IDs were valid.", "<br/>",
        length(genes_of_interest_found), " genes correctly mapped to UniProt IDs."))
    })
    
    output$numIsoformsMappedUP <- renderUI({
      if(input$no_isoforms){
        HTML(paste0("Alternative isoforms were excluded. ", nrow(ph_filt[ph_filt$gene %in% genes_of_interest,]), " proteins are shown. Isoform with earliest emergence shown."))
      }else{
        HTML(paste0("Alternative isoforms were included. ", nrow(ph[ph$Gene %in% genes_of_interest,]), " protein isoforms are shown. Any isoform emergence shown."))
      }
    })
    
    genes_of_interest
  })
  
  ##### disease_name #####
  disease_name <- reactive({
    if(!is.null(input$other_disease_of_interest) & input$runmode){
      input$other_disease_of_interest
    }else{
      input$disease_of_interest
    }
  })
  
  observeEvent(input$saverds, {
    if(is.na(input$other_disease_of_interest)){
      cache_name <- input$custom_file$datapath
    }else{
      cache_name <- input$other_disease_of_interest
    }
    
    MasterGeneLists[[cache_name]] = c(genes_of_interest_reactive(), rep(NA, nrow(MasterGeneLists) - length(genes_of_interest_reactive())))
    
    saveRDS(MasterGeneLists, file = "data/MasterGeneLists_Updated.rds")
    output$cache_message <- renderText({"Successfully cached data!"})
  })
  
  ##### phylotable #####
  output$phylo_table <- DT::renderDataTable({
    #input$run
    if(input$no_isoforms){
      ph <- ph_filt
    }
    
    if(input$show_uniprot_ids){
      ph <- ph %>%
        filter(gene %in% genes_of_interest_reactive()) %>%
        #filter(mrca_name %in% input$ps) %>%
        dplyr::select(Gene = gene, Phylostrata = mrca_name, `Isoform (UniProt)` = qseqid) %>%
        mutate(Phylostrata = factor(Phylostrata, levels = levels(ph$mrca_name)))
    }else{
      ph <- ph %>%
        filter(gene %in% genes_of_interest_reactive()) %>%
        #filter(mrca_name %in% input$ps) %>%
        group_by(Gene = gene, Phylostrata = mrca_name) %>%
        summarize(`Total Isoforms` = n()) %>%
        mutate(Phylostrata = factor(Phylostrata, levels = levels(ph$mrca_name)))
    }

    ph
  }, 
  server = T,
  extensions = 'Buttons', 
  options = list(
    paging = TRUE,
    searching = TRUE,
    fixedColumns = TRUE,
    scrollX = TRUE,
    scrollY = 400,
    autoWidth = TRUE,
    ordering = TRUE,
    lengthMenu = c(10, 25, 50, 100),
    dom = 'Bfrtip',
    buttons = c('copy', 'csv', 'excel', "pageLength")
  )
  )
  
  ##### plot of PS emergence ##### 
  output$ph_plot <- renderPlot({
    # Filter alt isoforms if specified
    if(input$no_isoforms){
      ph <- ph_filt
    }
    # Set up ph object
    genes_of_interest <- genes_of_interest_reactive()
    ph$is_of_interest <- ph$Gene %in% genes_of_interest
    ph$annot_of_interest <-
      ifelse(ph$is_of_interest,
             disease_name(),
             "All Human Genes")
    ph$Gene <-
      gene_uniprot_conv_df$Gene[match(ph$qseqid, gene_uniprot_conv_df$Entry)]

    # Reframe data
    ph_by_ps <- ph %>%
      dplyr::distinct() %>%
      dplyr::mutate(
        mrca_name = as.factor(mrca_name),
        annot_of_interest = as.factor(annot_of_interest)
      ) %>%
      dplyr::count(mrca_name = as.factor(mrca_name),
                   annot_of_interest,
                   .drop = F) %>%
      dplyr::mutate(gene_relevance = annot_of_interest) %>%
      dplyr::mutate(gene_relevance = factor(gene_relevance, levels = c(disease_name(), "All Human Genes")))
    
    # Generate plot
    ggplot(
      data = ph_by_ps,
      aes(
        x = mrca_name,
        y = log10(n),
        group = gene_relevance,
        fill = gene_relevance,
        label = round(n, digits = 2)
      )
    ) +
      geom_line(aes(color = gene_relevance), linewidth = 1.2) +
      geom_label(
        label.size = 0.05,
        alpha = 0.3,
        size = 2.5,
        vjust = -0.25,
        show.legend = F, 
        nudge_y = 0.1, 
        #fontface = "bold", #size = 6
      ) +
      labs(
        y = bquote(Number ~ of ~ novel ~ genes ~ log[10] ~ scaled),
        x = "Phylostrata",
        fill = "",
        color = "",
        title = paste0(
          "Emergence of <span style='color:#E41A1C;'><strong>",
          disease_name(),
          "</strong></span> vs <span style='color:#377EB8;'><strong>Overall Human Genes</strong></span><br>"
        )
        # title = paste0(
        #   "Normalized odds ratio of emergence of <br><span style='color:#E41A1C;'><strong>",
        #   disease_name(),
        #   "</strong></span> to overall human genes<br>"
        # )
      ) +
      # ggtitle(paste0(
      #   "Emergence of \"",
      #   disease_name(),
      #   "\" vs All Human Genes\n\n\n"
      # ), ) +
      scale_color_brewer(palette = "Set1") +
      scale_fill_brewer(palette = "Set1") +
      ylim(c(0, 5)) +
      ggpubr::labs_pubr() +
      my_theme +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), 
            panel.background = element_blank(),
            plot.title.position = "plot",
            plot.title = ggtext::element_markdown(size = 12, lineheight = 1, hjust = 0.5),
            axis.text = element_text(size = 10, face = "bold"),
            legend.position = "none")# +
      # ggpubr::theme_pubr(x.text.angle = 90,
      #                    margin = T,
      #                    base_size = 11) +
      #
    
  })
  
  #### ph plot click ####
  output$selected_ph <- DT::renderDataTable({
    #input$ph_plot_click
    if(is.null(input$ph_plot_click)){
      return(NULL)
    }
    click <- input$ph_plot_click
    if(is.null(click$x) | is.null(click$y)){
      return(NULL)
    }
    x <- max(round(click$x), 1)
    y <- click$y
    x <- levels(ph$mrca_name)[x]
    y <- 10^y
    
    if(input$no_isoforms){
      ph_filt %>%
        filter(mrca_name == x & gene %in% genes_of_interest_reactive()) %>%
        dplyr::select(Gene = gene, `UniProt ID` = qseqid, `Phylostrata` = mrca_name)
    }else{
      ph %>%
        filter(mrca_name == x & gene %in% genes_of_interest_reactive()) %>%
        dplyr::select(Gene = gene, `UniProt ID` = qseqid, `Phylostrata` = mrca_name)
    }
  }, server = T)
  
  ##### string plot #####
  output$string_plot <- renderPlot({
    if(input$no_isoforms){
      ph <- ph_filt
    }
    string_db <- STRINGdb::STRINGdb$new(species = 9606, version = "12", score_threshold = 200, input_directory = ".")
    ph$is_of_interest <- ph$Gene %in% genes_of_interest_reactive()
    ph$annot_of_interest <-
      ifelse(ph$is_of_interest,
             disease_name(),
             "All Human Genes")
    
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      paste0(hcl(h = hues, l = 65, c = 100)[1:n], "FF")
    }
    
    withProgress(message = "Pulling data from STRINGdb...", value = 0, {
      ph_mapped <- string_db$map(ph %>% filter(is_of_interest & mrca_name %in% input$ps) %>% as.data.frame(), "gene", removeUnmappedRows = T)
      #payload_id <- string_db$post_payload( ph_mapped$STRING_id, colors=gg_color_hue(length(unique(ph_mapped$mrca_name))))
      incProgress(amount = 0.5, message = "Generating plot...")
      string_db$plot_network(ph_mapped$STRING_id, required_score = 800)
      incProgress(amount = 0.5, message = "Done!")
    })
  })
  
  ##### normalized evo rate  #####
  output$rate_plot <- renderPlot({
    if(input$no_isoforms){
      ph <- ph_filt
    }
    ph$is_of_interest <- ph$Gene %in% genes_of_interest_reactive()
    ph$annot_of_interest <-
      ifelse(ph$is_of_interest,
             disease_name(),
             "All Human Genes")

    ph_by_ps <- ph %>%
      dplyr::distinct() %>%
      dplyr::mutate(
        mrca_name = as.factor(mrca_name),
        annot_of_interest = relevel(as.factor(annot_of_interest), ref = "All Human Genes")
      ) %>%
      dplyr::count(mrca_name = as.factor(mrca_name),
                   annot_of_interest,
                   .drop = F) %>%
      dplyr::mutate(gene_relevance = annot_of_interest)
    
    
    fixed_ratio <- ph_by_ps %>%
      group_by(gene_relevance) %>%
      summarize(total_in_group = sum(n)) %>%
      pull(total_in_group)
    fixed_ratio <- min(fixed_ratio) / sum(fixed_ratio)
    
    ratio_plot_data <- ph_by_ps %>%
      pivot_wider(
        id_cols = c(mrca_name),
        names_from = gene_relevance,
        values_from = n
      ) %>%
      mutate(`Other All Human Genes` = sum(`All Human Genes`) - `All Human Genes`) %>%
      mutate(`Other Disease` = sum(.[[3]]) - .[[3]]) %>%
      mutate(`Odds Ratio` = (.[[3]] / `Other Disease`) / .[[2]] / (`Other All Human Genes`)) %>%
      mutate(Expected = fixed_ratio * `All Human Genes`) %>%
      mutate(`ChiSq` = (.[[3]] - Expected)^2 / Expected) %>%
      mutate(pchisq = pchisq(ChiSq, 1, lower.tail = F)) %>%
      mutate(Ratio = (.[[3]] / .[[2]]) / fixed_ratio) %>%
      filter(!is.na(Ratio))
    
    ggplot(
      data = ratio_plot_data,
      aes(x = mrca_name, y = Ratio)
    ) +
      geom_hline(yintercept = 1, color = "grey80") +
      geom_point() +
      #geom_ribbon(aes(ymin = 1, ymax = pmax(Ratio, 1), x = mrca_name), fill = "#E41A1C", group = "a") +
      #ggpubr::theme_pubr(x.text.angle = 90) +
      geom_line(group = "a", linewidth = 1.25) +
      geom_text(data = ratio_plot_data %>% filter(Ratio >= 1),
                aes(label = ifelse(pchisq < 0.05, "*", "")), vjust = 1, nudge_y = 1, size = 5, color = "red") +
      labs(
        x = "Phylostrata",
        y = paste0(
          "Odds Ratio of Evolutionary Rate \n(Disease Genes : Total Human Genome)"
        ),
        title = paste0(
          "Normalized odds ratio of emergence of <br><span style='color:#E41A1C;'><strong>",
          disease_name(),
          "</strong></span> to overall human genes<br>"
        )
      ) +
      ggpubr::labs_pubr() +
      my_theme +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            plot.title.position = "plot",
            plot.title = ggtext::element_markdown(size = 12, lineheight = 1, hjust = 0.5))
  })
  
  ##### GO reactive data ####
  GO_df_reactive <- reactive({
    input$run
    
    if(isolate(input$no_isoforms)){
      ph <- ph_filt
    }
    ph$is_of_interest <- ph$Gene %in% genes_of_interest_reactive()
    ph$annot_of_interest <-
      ifelse(ph$is_of_interest,
             disease_name(),
             "All Human Genes")

    ph_of_interest <- ph %>%
      filter(is_of_interest == TRUE)
    
    gene_list_by_strata <-
      lapply(unique(ph_of_interest$ps), function(x) {
        ph_of_interest$gene[ph_of_interest$ps == x]
      })
    names(gene_list_by_strata) <- unique(ph_of_interest$ps)
    
    # This takes a while, but these are the best results...
    annFUN.org.symbol <-
      topGO::annFUN.org(whichOnto = "BP",
                        mapping = "org.Hs.eg.db",
                        ID = "symbol")
    
    GO_df_by_strata <-
      lapply(names(gene_list_by_strata), function(i) {
        gene_set_to_enrich <- unique(gene_list_by_strata[[i]])
        gene_set_to_enrich <-
          gene_set_to_enrich[!is.na(gene_set_to_enrich)]
        if (length(gene_set_to_enrich) < 1) {
          return(NA)
        }
        
        allGenes <- unique(unlist(annFUN.org.symbol))
        geneList <-
          as.factor(as.integer(allGenes %in% gene_set_to_enrich))
        names(geneList) <- allGenes
        geneList_numeric <- as.numeric(geneList[geneList == 1])
        names(geneList_numeric) <- names(geneList)[geneList == 1]
        
        pathways <- annFUN.org.symbol
        
        # NEW WAY USING FGSEA -
        allres <- fgsea(pathways = pathways,
                        stats = geneList_numeric,
                        minSize = 1)
        
        # OLD WAY USING TOPGO - 
        # GOdata <- new(
        #   "topGOdata",
        #   ontology = "BP",
        #   allGenes = geneList,
        #   annot = annFUN.org,
        #   mapping = "org.Hs.eg.db",
        #   ID = "symbol"
        # )
        # GOgraph <- graph(GOdata)
        # gs <- geneScore(GOdata, whichGenes = gene_set_to_enrich)
        
        # resultFisher <-
        #   runTest(GOdata, algorithm = "classic", statistic = "fisher")
        # resultKS <- runTest(GOdata, algorithm = "elim", statistic = "ks")
        
        #allres <- GenTable(GOdata, classic = resultFisher, numChar = 1000)
        ## Add progress:
        #incProgress(1/length(gene_list_by_strata), detail = paste0("Computing enrichment for ", i))
        
        return(allres)
        
        # return(
        #   res %>%
        #     dplyr::mutate(
        #       Pathway_Name = gsub(
        #         pattern = "GO(BP|MF|CC)",
        #         x = gsub("_", x = pathway, replacement = " "),
        #         replacement = ""
        #       )
        #     ) %>%
        #     dplyr::mutate(Pathway_Name = str_to_title(Pathway_Name)) %>%
        #     dplyr::slice_max(abs(NES), n = 30)
        # )
      }) 
    
    
    names(GO_df_by_strata) <-
      ph$mrca_name[match(names(gene_list_by_strata), ph$ps)]
    
    GO_res_df <-
      bind_rows(lapply(GO_df_by_strata, as.data.frame), .id = "PS") %>%
      mutate(PS = factor(PS, levels = levels(ph$mrca_name))) %>%
      mutate(Term = AnnotationDbi::Term(pathway)) %>%
      filter(!is.na(pval))
    
    GO_res_df
  })
  
  ##### Dup GO term plot #####
  #output$duplicated_go_plot <- plotly::renderPlotly({
    # ph$is_of_interest <- ph$qseqid %in% up_of_interest_reactive()
    # ph$annot_of_interest <-
    #   ifelse(ph$is_of_interest,
    #          input$disease_of_interest,
    #          "All Human Genes")
    # ph$Gene <-
    #   gene_uniprot_conv_df$To[match(ph$qseqid, gene_uniprot_conv_df$From)]
    # 
    # ph_of_interest <- ph %>%
    #   filter(is_of_interest == TRUE)
    # 
    # gene_list_by_strata <-
    #   lapply(unique(ph_of_interest$ps), function(x) {
    #     ph_of_interest$gene[ph_of_interest$ps == x]
    #   })
    # names(gene_list_by_strata) <- unique(ph_of_interest$ps)
    # 
    # annFUN.org.symbol <-
    #   topGO::annFUN.org(whichOnto = "BP",
    #                     mapping = "org.Hs.eg.db",
    #                     ID = "symbol")
    # 
    # GO_df_by_strata <-
    #   lapply(names(gene_list_by_strata), function(i) {
    #     gene_set_to_enrich <- unique(gene_list_by_strata[[i]])
    #     gene_set_to_enrich <-
    #       gene_set_to_enrich[!is.na(gene_set_to_enrich)]
    #     if (length(gene_set_to_enrich) < 2) {
    #       return(NA)
    #     }
    #     
    #     allGenes <- unique(unlist(annFUN.org.symbol))
    #     geneList <-
    #       as.factor(as.integer(allGenes %in% gene_set_to_enrich))
    #     names(geneList) <- allGenes
    #     geneList_numeric <- as.numeric(geneList[geneList == 1])
    #     names(geneList_numeric) <- names(geneList)[geneList == 1]
    #     
    #     res <- fgsea(pathways = pathways,
    #                  stats = geneList_numeric,
    #                  minSize = 0)
    #     
    #     # GOdata <- new(
    #     #   "topGOdata",
    #     #   ontology = "BP",
    #     #   allGenes = geneList,
    #     #   annot = annFUN.org,
    #     #   mapping = "org.Hs.eg.db",
    #     #   ID = "symbol"
    #     # )
    #     # GOgraph <- graph(GOdata)
    #     # gs <- geneScore(GOdata, whichGenes = gene_set_to_enrich)
    #     
    #     # resultFisher <-
    #     #   runTest(GOdata, algorithm = "classic", statistic = "fisher")
    #     # resultKS <- runTest(GOdata, algorithm = "elim", statistic = "ks")
    #     
    #     # allres <- GenTable(GOdata, classic = resultFisher, numChar = 1000)
    #     
    #     # return(allres)
    #     
    #     return(
    #       res %>%
    #         dplyr::mutate(
    #           Pathway_Name = gsub(
    #             pattern = "GO(BP|MF|CC)",
    #             x = gsub("_", x = pathway, replacement = " "),
    #             replacement = ""
    #           )
    #         ) %>%
    #         dplyr::mutate(Pathway_Name = str_to_title(Pathway_Name)) %>%
    #         dplyr::slice_max(abs(NES), n = 30)
    #     )
    #   })
    
  #   GO_res_df <- GO_df_reactive()
  #   
  #   if (nrow(GO_res_df) == 0) {
  #     plot_out <- ggplot()
  #   } else{
  #     GO_res_df_dup <- GO_res_df %>%
  #       dplyr::filter(base::duplicated(Term) |
  #                       base::duplicated(Term, fromLast = T)) %>%
  #       dplyr::filter(!is.na(Term)) %>%
  #       dplyr::filter(pval < 0.1)
  #     
  #     # IF too many terms
  #     if(nrow(GO_res_df_dup) > 20){
  #       # Get the terms that appear the most number of times
  #       term_counts <- table(GO_res_df_dup$Term)
  #       top20_terms <- names(term_counts)[order(term_counts, decreasing = T)[1:20]]
  #       GO_res_df_dup <- GO_res_df_dup %>% filter(Term %in% top20_terms)
  #     }
  #     
  #     plot_out <- ggplot(GO_res_df_dup,
  #                        aes(x = PS, y = Term, group = Term)) +
  #       geom_line() +
  #       geom_point(aes(color = -log10(as.numeric(pval)))) +
  #       #scale_x_continuous(n.breaks = 27) +
  #       ggpubr::theme_pubr(x.text.angle = 90, base_size = 8, legend = "none") +
  #       labs(
  #         x = "Phylostrata",
  #         y = "GO Term",
  #         title = paste0(
  #           "BPs with multiple evolutionary origins related to ",
  #           disease_name()
  #         ),
  #         color = "Hypergeometric Test -log10 P-value"
  #       )
  #     
  #     if (nrow(GO_res_df_dup) == 0) {
  #       plot_out <- plot_out +
  #         annotate(
  #           geom = "text",
  #           x = 1,
  #           y = 1,
  #           label = "No repeated GO terms over phylostrata."
  #         )
  #     }
  #   }
  #   
  #   plotly::ggplotly(p = plot_out)
  # })
  
  ##### Plotly GO #####
  # output$GOplot2 <- plotly::renderPlotly({
  #   GO_res_df <- GO_df_reactive()
  #   
  #   GO_res_df_top <- GO_res_df %>%
  #     group_by(PS) %>%
  #     arrange(desc(padj)) %>%
  #     #slice_min(n = input$maxGoTerms, order_by = classic) %>%
  #     ungroup()
  #   
  #   plot2 <- ggplot(GO_res_df_top, aes(
  #     x = PS,
  #     y = -log10(as.numeric(pval)),
  #     label = `Term`,
  #     color = PS
  #   )) +
  #     # geom_hline(
  #     #   yintercept = -log10(0.05),
  #     #   linetype = 2,
  #     #   color = "grey70"
  #     # ) +
  #     geom_point(show.legend = F, position = position_jitter(0.1)) +
  #     geom_label(
  #       data = GO_res_df_top,
  #       aes(label = `Term`, y = -log(as.numeric(pval))),
  #       show.legend = F,
  #       size = 3,
  #       hjust = 0.5,
  #       vjust = 1
  #     ) +
  #     scale_x_discrete(limits = levels(GO_res_df$PS)) +
  #     ggpubr::theme_pubr(x.text.angle = 90, legend = "none") +
  #     labs(
  #       x = "Phylostrata (PS)",
  #       y = "-log P-value (Fisher test)",
  #       title = "Top enriched GO terms by phylostrata",
  #       subtitle = "Top 3 in each PS are labeled"
  #     )
  # 
  #   plotly::ggplotly(plot2, tooltip = "label")
  # })
  
    ##### Static GO ####
  output$GOplot_static <- renderPlot({
    GO_res_df <- GO_df_reactive()
    
    cappval <- isolate(input$cappval)
    if(cappval < 100){
      GO_res_df$pval[-log10(GO_res_df$pval) > cappval] <- 10^(-1*input$cappval)
    }
    
    min_lim <- isolate(input$pval_range[1])
    min_lim <- ifelse(min_lim < 10^(-1*cappval), 10^(-1*cappval), min_lim)
    
    max_lim <- isolate(input$pval_range[2])
    max_lim <- ifelse(max_lim == min_lim, 1, max_lim)
    
    GO_res_df <- GO_res_df %>%
      dplyr::filter(pval >= min_lim & pval <= max_lim)
    
    
    GO_res_df_top <- GO_res_df %>%
      group_by(PS) %>%
      arrange(desc(pval)) %>%
      slice_min(n = isolate(input$labeltop), order_by = pval, with_ties = F) %>%
      ungroup()
    
    palette <- scales::hue_pal()(length(unique(ph$mrca_name)))
    names(palette) = levels(ph$mrca_name)
    
    ggplot(GO_res_df, aes(
      x = PS,
      y = as.numeric(pval),
      label = `Term`,
      color = PS
    )) +
      geom_point(show.legend = F, position = position_jitter(0.1), size = 4) +
      ggrepel::geom_label_repel(
        data = GO_res_df_top,
        aes(label = `Term`, y = as.numeric(pval)),
        show.legend = F,
        size = 3,
        hjust = 0.5,
        vjust = 1, max.overlaps = 30
      ) +
      scale_color_manual(values = palette) +
      scale_x_discrete(limits = levels(GO_res_df$PS)) +
      scale_y_continuous(transform = scales::new_transform(name = "-log10", 
                                                           transform = function(x) -log10(x), 
                                                           inverse = function(x) 10^(-1*x))) +
      # scale_y_continuous(transform = "log10", 
      #                    breaks = c(0.999, 0.99, 0.95, 0.5, 0.01), 
      #                    labels = c("0.001", "0.01", "0.05", "0.5", ">0.99")) +
      # ggpubr::theme_pubr(x.text.angle = 90, legend = "none", 
      #                    base_family = "arial", base_size = 16, border = T) +
      my_theme +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(
        x = "Phylostrata (PS)",
        y = "P-value (Fisher test)",
        title = paste0("Top enriched GO terms by phylostrata for ", disease_name())
      ) +
      ggpubr::labs_pubr()
    #plotly::ggplotly(plot2, tooltip = "label")
  })
  
  #### GO plot click ####
  output$selected_go <- renderPlot({
    #input$ph_plot_click
    if(is.null(input$go_plot_click)){
      return(NULL)
    }
    click <- input$go_plot_click
    if(is.null(click$x) | is.null(click$y)){
      return(NULL)
    }
    x <- max(round(click$x), 1)
    y <- click$y
    x <- levels(ph$mrca_name)[x]
    
    palette <- scales::hue_pal()(length(unique(ph$mrca_name)))
    names(palette) = levels(ph$mrca_name)
    
    GO_res_df <- GO_df_reactive()
    GO_res_df <- GO_res_df %>%
      filter(PS == x) %>%
      dplyr::select(Term = Term, `P-value` = pval, `Normalized Enrichment Score` = NES, `Adjusted P-value` = padj)
    
    ggplot(data = GO_res_df, aes(x = `Normalized Enrichment Score`, y = `Adjusted P-value`)) +
      geom_point(color = palette[x]) +
      scale_color_manual(values = palette[x]) +
      ggrepel::geom_label_repel(
        data = GO_res_df %>% slice_min(`Adjusted P-value`, n = input$labeltop, with_ties = F),
        aes(label = Term), max.overlaps = 30) +
      #scale_x_continuous(limits = c(-2, 2)) +
      scale_y_continuous(transform = scales::new_transform(name = "-log10", 
                                                           transform = function(x) -log10(x), 
                                                           inverse = function(x) 10^(-1*x))) +
      my_theme +
      labs(
        x = "Normalized Enrichment Score",
        y = "-log10 Adjusted P-value",
        title = paste0("Top enriched GO terms for ", disease_name(), " \nin PS ", x),
      ) +
      ggpubr::labs_pubr()
    
  })
  
  ##### string plot #####
  output$string_plot <- renderPlot({
    if(input$no_isoforms){
      ph <- ph_filt
    }
    string_db <- STRINGdb::STRINGdb$new(species = 9606, version = "12", score_threshold = 200, input_directory = ".")
    ph$is_of_interest <- ph$Gene %in% genes_of_interest_reactive()
    ph$annot_of_interest <-
      ifelse(ph$is_of_interest,
             disease_name(),
             "All Human Genes")
    
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      paste0(hcl(h = hues, l = 65, c = 100)[1:n], "FF")
    }
    
    withProgress(message = "Pulling data from STRINGdb...", value = 0, {
      ph_mapped <- string_db$map(ph %>% filter(is_of_interest & mrca_name %in% input$ps) %>% as.data.frame(), "gene", removeUnmappedRows = T)
      #payload_id <- string_db$post_payload( ph_mapped$STRING_id, colors=gg_color_hue(length(unique(ph_mapped$mrca_name))))
      incProgress(amount = 0.5, message = "Generating plot...")
      string_db$plot_network(ph_mapped$STRING_id, required_score = 800)
      incProgress(amount = 0.5, message = "Done!")
    })
  })
  
}

##### RUN #####
shinyApp(ui = ui, server = server)
