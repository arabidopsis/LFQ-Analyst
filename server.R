# Define server logic to read selected file ----
library("shiny", quietly = TRUE)
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 100 * 1024^2) ## Set maximum upload size to 100MB

  #  Show elements on clicking Start analysis button
  observeEvent(input$analyze, {
    if (input$analyze == 0) {
      return()
    }
    shinyjs::hide("quickstart_info")
    shinyjs::show("downloadbox")
  })

  observeEvent(input$analyze, {
    if (input$analyze == 0) {
      return()
    }
    shinyjs::show("results_tab")
  })

  observeEvent(input$analyze, {
    if (input$analyze == 0) {
      return()
    }
    shinyjs::show("qc_tab")
  })

  observeEvent(input$analyze, {
    if (input$analyze == 0) {
      return()
    }
    shinyjs::show("enrichment_tab")
  })

  ## Shinyalert
  observeEvent(input$analyze, {
    if (input$analyze == 0) {
      return()
    }

    shinyalert::shinyalert("In Progress!", "Data analysis has started, wait until table and plots
                appear on the screen",
      type = "info",
      closeOnClickOutside = TRUE,
      closeOnEsc = TRUE,
      timer = 10000
    ) # timer in miliseconds (10 sec)
  })


  #### ======= Render Functions

  output$volcano_cntrst <- renderUI({
    if (!is.null(comparisons())) {
      df <- SummarizedExperiment::rowData(dep())
      cols <- grep("_significant$", colnames(df))
      selectizeInput("volcano_cntrst",
        "Comparison",
        choices = gsub("_significant", "", colnames(df)[cols])
      )
    }
  })

  ## comparisons
  output$contrast <- renderUI({
    if (!is.null(comparisons())) {
      df <- SummarizedExperiment::rowData(dep())
      cols <- grep("_significant$", colnames(df))
      selectizeInput("contrast",
        "Comparison",
        choices = gsub("_significant", "", colnames(df)[cols])
      )
    }
  })

  output$contrast_1 <- renderUI({
    if (!is.null(comparisons())) {
      df <- SummarizedExperiment::rowData(dep())
      cols <- grep("_significant$", colnames(df))
      selectizeInput("contrast",
        "Comparison",
        choices = gsub("_significant", "", colnames(df)[cols])
      )
    }
  })

  output$downloadTable <- renderUI({
    if (!is.null(dep())) {
      selectizeInput(
        "dataset",
        "Choose a dataset to save",
        c(
          "Results", "Original_matrix",
          "Imputed_matrix",
          "Full_dataset"
        )
      )
    }
  })

  output$downloadButton <- renderUI({
    if (!is.null(dep())) {
      downloadButton("downloadData", "Save")
    }
  })

  output$downloadZip <- renderUI({
    if (!is.null(dep())) {
      downloadButton("downloadZip1", "Download result plots")
    }
  })
  output$downloadreport <- renderUI({
    if (!is.null(dep())) {
      downloadButton("downloadReport", "Download Report")
    }
  })

  output$downloadPlots <- renderUI({
    if (!is.null(dep())) {
      downloadButton("downloadPlots1", "Download Plots")
    }
  })


  ## make reactive elements
  maxquant_data_input <- reactive({
    NULL
  })
  exp_design_input <- reactive({
    NULL
  })


  maxquant_data_input <- eventReactive(input$analyze, {
    inFile <- input$file1
    if (is.null(inFile)) {
      return(NULL)
    }
    temp_data <- read.table(inFile$datapath,
      header = TRUE,
      fill = TRUE, # to fill any missing data
      sep = "\t",
      quote = ""
    )
    validate(maxquant_input_test(temp_data))
    return(temp_data)
  })

  # observeEvent(input$analyze,{
  #   exp_design<-reactive({
  #     inFile<-input$file2
  #     if (is.null(inFile))
  #       return(NULL)
  #     temp_df<-read.table(inFile$datapath,
  #                         header = TRUE,
  #                         sep="\t",
  #                         stringsAsFactors = FALSE)
  #     exp_design_test(temp_df)
  #     temp_df$label<-as.character(temp_df$label)
  #     return(temp_df)
  #   })
  # })
  exp_design_input <- eventReactive(input$analyze, {
    inFile <- input$file2
    if (is.null(inFile)) {
      return(NULL)
    }
    temp_df <- read.table(inFile$datapath,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE
    )
    exp_design_test(temp_df)
    temp_df$label <- as.character(temp_df$label)
    temp_df$condition <- trimws(temp_df$condition, which = "left")
    temp_df$condition <- temp_df$condition %>% gsub("[^[:alnum:]|_]+", "_", .) # auto fix special characters
    return(temp_df)
  })



  ### Reactive components
  processed_data <- reactive({
    ## check which dataset
    if (!is.null(maxquant_data_input())) {
      maxquant_data <- reactive({
        maxquant_data_input()
      })
    }

    if (!is.null(exp_design_input())) {
      exp_design <- reactive({
        exp_design_input()
      })
    }


    message(exp_design())
    if (any(grepl("+", maxquant_data()$Reverse))) {
      filtered_data <- dplyr::filter(maxquant_data(), Reverse != "+")
    } else {
      filtered_data <- maxquant_data()
    }
    if (any(grepl("+", filtered_data$Potential.contaminant))) {
      filtered_data <- dplyr::filter(filtered_data, Potential.contaminant != "+")
    }
    if (any(grepl("+", filtered_data$Only.identified.by.site))) {
      filtered_data <- dplyr::filter(filtered_data, Only.identified.by.site != "+")
    }
    if (input$single_peptide == TRUE) {
      filtered_data <- filtered_data
    } else {
      filtered_data <- dplyr::filter(filtered_data, as.numeric(Razor...unique.peptides) >= 2)
    }

    filtered_data <- ids_test(filtered_data)

    data_unique <- DEP::make_unique(filtered_data, "Gene.names", "Protein.IDs", delim = ";")
    lfq_columns <- grep("LFQ.", colnames(data_unique))

    # ensure all intensity columns are numeric type
    data_unique[, lfq_columns] <- sapply(data_unique[, lfq_columns], as.numeric)
    ## Check for matching columns in maxquant and experiment design file
    test_match_lfq_column_design(data_unique, lfq_columns, exp_design())
    data_se <- DEP:::make_se(data_unique, lfq_columns, exp_design())

    # # Check number of replicates
    # if(max(exp_design()$replicate)<3){
    #   threshold<-0
    # } else  if(max(exp_design()$replicate)==3){
    #   threshold<-1
    # } else if(max(exp_design()$replicate)<6 ){
    #   threshold<-2
    # } else if (max(exp_design()$replicate)>=6){
    #   threshold<-trunc(max(exp_design()$replicate)/2)
    # }
    #
    #
    # filter_missval(data_se,thr = threshold)
    # filter_missval(data_se,thr = threshold)
    exp_df <- exp_design() %>% dplyr::count(condition)
    exp_df <- exp_df %>% dplyr::mutate(thr = lapply(exp_df$n, threshold_detect)) # function:threshold_detect
    condition_list <- exp_df$condition

    data_filtered <- filter_missval_new(data_se, condition_list, exp_df)
    return(data_filtered)
  })

  unimputed_table <- reactive({
    temp <- SummarizedExperiment::assay(processed_data())
    temp1 <- 2^(temp)
    colnames(temp1) <- paste(colnames(temp1), "original_intensity", sep = "_")
    temp1 <- cbind(ProteinID = rownames(temp1), temp1)
    # temp1$ProteinID<-rownames(temp1)
    return(as.data.frame(temp1))
  })

  normalised_data <- reactive({
    DEP::normalize_vsn(processed_data())
  })

  imputed_data <- reactive({
    DEP::impute(processed_data(), input$imputation)
  })

  imputed_table <- reactive({
    temp <- SummarizedExperiment::assay(imputed_data())
    # tibble::rownames_to_column(temp,var = "ProteinID")
    temp1 <- 2^(temp)
    colnames(temp1) <- paste(colnames(temp1), "imputed_intensity", sep = "_")
    temp1 <- cbind(ProteinID = rownames(temp1), temp1) # temp1$ProteinID<-rownames(temp1)
    return(as.data.frame(temp1))
  })

  diff_all <- reactive({
    DEP::test_diff(imputed_data(), type = "all")
  })

  dep <- reactive({
    if (input$fdr_correction == "BH") {
      diff_all <- test_limma(imputed_data(), type = "all", paired = input$paired)
      DEP::add_rejections(diff_all, alpha = input$p, lfc = input$lfc)
    } else {
      diff_all <- test_diff(imputed_data(), type = "all")
      DEP::add_rejections(diff_all, alpha = input$p, lfc = input$lfc)
    }
  })

  comparisons <- reactive({
    temp <- capture.output(DEP::test_diff(imputed_data(), type = "all"), type = "message")
    gsub(".*: ", "", temp)
    ## Split conditions into character vector
    unlist(strsplit(temp, ","))
    ## Remove leading and trailing spaces
    trimws(temp)
  })


  ## Results plot inputs

  ## PCA Plot
  pca_input <- eventReactive(input$analyze, {
    if (input$analyze == 0) {
      return()
    }
    if (num_total() <= 500) {
      if (length(levels(as.factor(SummarizedExperiment::colData(dep())$replicate))) <= 6) {
        pca_plot <- DEP::plot_pca(dep(), n = num_total(), point_size = 4)
        pca_plot <- pca_plot + labs(title = "PCA Plot")
        return(pca_plot)
      } else {
        pca_plot <- DEP::plot_pca(dep(), n = num_total(), point_size = 4, indicate = "condition")
        pca_plot <- pca_plot + labs(title = "PCA Plot")
        return(pca_plot)
      }
    } else {
      if (length(levels(as.factor(SummarizedExperiment::colData(dep())$replicate))) <= 6) {
        pca_plot <- DEP::plot_pca(dep(), point_size = 4)
        pca_plot <- pca_plot + labs(title = "PCA Plot")
        return(pca_plot)
      } else {
        # pca_label<-SummarizedExperiment::colData(dep())$replicate
        pca_plot <- DEP::plot_pca(dep(), point_size = 4, indicate = "condition")
        # pca_plot<-pca_plot + geom_point()
        pca_plot <- pca_plot + ggrepel::geom_text_repel(aes(label = factor(rowname)),
          size = 4,
          box.padding = unit(0.1, "lines"),
          point.padding = unit(0.1, "lines"),
          segment.size = 0.5
        )
        pca_plot <- pca_plot + labs(title = "PCA Plot")

        #        pca_plot<-DEP::plot_pca(dep(), point_size = 4, indicate = "condition")
        #         pca_plot + ggrepel::geom_text_repel(aes(label=SummarizedExperiment::colData(dep())$replicate),
        #                                        size = 5,
        #                                           box.padding = unit(0.1, 'lines'),
        #                                          point.padding = unit(0.1, 'lines'),
        #                                         segment.size = 0.5)
        return(pca_plot)
      }
    }
  })

  ### Heatmap Differentially expressed proteins
  heatmap_cluster <- eventReactive(input$analyze, {
    if (input$analyze == 0) {
      return()
    }
    heatmap_list <- get_cluster_heatmap(dep(),
      type = "centered", kmeans = TRUE,
      k = input$k_number, col_limit = 6,
      indicate = "condition"
    )
    return(heatmap_list)
  })

  heatmap_input <- reactive({
    heatmap_list <- heatmap_cluster()
    heatmap_list[[1]]
  })


  # heatmap_input<-eventReactive(input$analyze ,{
  #   if(input$analyze==0 ){
  #     return()
  #   }
  #
  #   heatmap_list <-get_cluster_heatmap(dep(),
  #                       type="centered",kmeans = TRUE,
  #                       k=input$k_number, col_limit = 6,
  #                       indicate = "condition"
  #                       )
  #   heatmap_list[[1]]
  # })

  ### Volcano Plot
  volcano_input <- reactive({
    if (!is.null(input$volcano_cntrst)) {
      plot_volcano_new(
        dep(),
        input$volcano_cntrst,
        input$fontsize,
        input$check_names,
        input$p_adj
      )
    }
  })

  volcano_df <- reactive({
    if (!is.null(input$volcano_cntrst)) {
      get_volcano_df(
        dep(),
        input$volcano_cntrst,
        input$p_adj
      )
    }
  })


  volcano_input_selected <- reactive({
    if (!is.null(input$volcano_cntrst)) {
      if (!is.null(input$contents_rows_selected)) {
        proteins_selected <- data_result()[c(input$contents_rows_selected), ] ## get all rows selected
      } else if (!is.null(input$protein_brush)) {
        proteins_selected <- data_result()[data_result()[["Gene Name"]] %in% protein_name_brush(), ]
      }

      ## convert contrast to x and padj to y
      diff_proteins <- grep(
        paste("^", input$volcano_cntrst, "_log2", sep = ""),
        colnames(proteins_selected)
      )
      if (input$p_adj == "FALSE") {
        padj_proteins <- grep(
          paste("^", input$volcano_cntrst, "_p.val", sep = ""),
          colnames(proteins_selected)
        )
      } else {
        padj_proteins <- grep(
          paste("^", input$volcano_cntrst, "_p.adj", sep = ""),
          colnames(proteins_selected)
        )
      }

      df_protein <- data.frame(
        x = proteins_selected[, diff_proteins],
        y = -log10(as.numeric(proteins_selected[, padj_proteins])), # )#,
        name = proteins_selected$`Gene Name`
      )
      # print(df_protein)
      p <- plot_volcano_new(
        dep(),
        input$volcano_cntrst,
        input$fontsize,
        input$check_names,
        input$p_adj
      )

      p + geom_point(data = df_protein, aes(x, y), color = "maroon", size = 3) +
        ggrepel::geom_text_repel(
          data = df_protein,
          aes(x, y, label = name),
          size = 4,
          box.padding = unit(0.1, "lines"),
          point.padding = unit(0.1, "lines"),
          segment.size = 0.5
        ) ## use the dataframe to plot points
    }
  })

  protein_input <- reactive({
    # protein_selected  <- data_result()[input$contents_rows_selected,1]
    df <- data_result()
    if (!is.null(input$protein_brush)) {
      df <- df[df[["Gene Name"]] %in% protein_name_brush(), ]
    }

    if (!is.null(input$protein_click)) {
      df <- df[df[["Gene Name"]] %in% protein_name_click(), ]
    }
    protein_selected <- df[input$contents_rows_selected, 1]

    if (length(levels(as.factor(SummarizedExperiment::colData(dep())$replicate))) <= 8) {
      plot_protein(dep(), protein_selected, input$type)
    } else {
      protein_plot <- plot_protein(dep(), protein_selected, input$type)
      protein_plot + scale_color_brewer(palette = "Paired")
    }
  })


  ## QC Inputs
  norm_input <- reactive({
    DEP::plot_normalization(
      processed_data(),
      normalised_data()
    )
  })

  missval_input <- reactive({
    DEP::plot_missval(processed_data())
  })

  detect_input <- reactive({
    DEP::plot_detect(processed_data())
  })

  imputation_input <- reactive({
    DEP::plot_imputation(
      normalised_data(),
      diff_all()
    )
  })

  p_hist_input <- reactive({
    plot_p_hist(dep())
  })

  numbers_input <- reactive({
    DEP::plot_numbers(normalised_data())
  })

  coverage_input <- reactive({
    DEP::plot_coverage(normalised_data())
  })

  correlation_input <- reactive({
    DEP::plot_cor(dep(), significant = FALSE)
  })

  cvs_input <- reactive({
    plot_cvs(dep())
  })

  num_total <- reactive({
    dep() %>%
      nrow()
  })

  ## Enrichment inputs

  go_input <- eventReactive(input$go_analysis, {
    withProgress(
      message = "Gene ontology enrichment is in progress",
      detail = "Please wait for a while",
      value = 0,
      {
        for (i in 1:15) {
          incProgress(1 / 15)
          Sys.sleep(0.25)
        }
      }
    )

    if (!is.null(input$contrast)) {
      enrichment_output_test(dep(), input$go_database)
      go_results <- test_gsea_mod(dep(), databases = input$go_database, contrasts = TRUE)
      null_enrichment_test(go_results)
      plot_go <- plot_enrichment(go_results,
        number = 5, alpha = 0.05, contrasts = input$contrast,
        databases = input$go_database, nrow = 2, term_size = 8
      ) + aes(stringr::str_wrap(Term, 60)) +
        xlab(NULL)
      go_list <- list("go_result" = go_results, "plot_go" = plot_go)
      return(go_list)
    }
  })

  pathway_input <- eventReactive(input$pathway_analysis, {
    progress_indicator("Pathway Analysis is running....")
    enrichment_output_test(dep(), input$pathway_database)
    pathway_results <- test_gsea_mod(dep(), databases = input$pathway_database, contrasts = TRUE)
    null_enrichment_test(pathway_results)
    plot_pathway <- plot_enrichment(pathway_results,
      number = 5, alpha = 0.05, contrasts = input$contrast_1,
      databases = input$pathway_database, nrow = 3, term_size = 8
    ) + aes(stringr::str_wrap(Term, 30)) +
      xlab(NULL)
    pathway_list <- list("pa_result" = pathway_results, "plot_pa" = plot_pathway)
    return(pathway_list)
  })


  #### Interactive UI
  output$significantBox <- shinydashboard::renderInfoBox({
    num_total <- dep() %>%
      nrow()
    num_signif <- dep() %>%
      .[SummarizedExperiment::rowData(.)$significant, ] %>%
      nrow()
    frac <- num_signif / num_total

    info_box <- shinydashboard::infoBox("Significant proteins",
      paste0(
        num_signif,
        " out of ",
        num_total
      ),
      paste0(
        signif(frac * 100, digits = 3),
        "% of proteins differentially expressed across all conditions"
      ),
      icon = icon("stats", lib = "glyphicon"),
      color = "olive",
      # fill = TRUE,
      width = 4
    )

    return(info_box)
  })

  ##### Get results dataframe from Summarizedexperiment object
  data_result <- reactive({
    get_results_proteins(dep())
    # get_results(dep())
  })


  #### Data table
  output$contents <- DT::renderDataTable(
    {
      df <- data_result()
      return(df)
    },
    options = list(
      scrollX = TRUE,
      autoWidth = TRUE,
      columnDefs = list(list(width = "400px", targets = c(-1)))
    )
  )

  ## Deselect all rows button
  proxy <- DT::dataTableProxy("contents")

  observeEvent(input$clear, {
    proxy %>% DT::selectRows(NULL)
  })

  observeEvent(input$original, {
    output$contents <- DT::renderDataTable(
      {
        df <- data_result()
        return(df)
      },
      options = list(
        scrollX = TRUE,
        autoWidth = TRUE,
        columnDefs = list(list(width = "400px", targets = c(-1)))
      )
    )
  })

  protein_name_brush <- reactive({
    # protein_tmp<-nearPoints(volcano_df(), input$protein_click, maxpoints = 1)
    protein_tmp <- brushedPoints(volcano_df(), input$protein_brush,
      xvar = "diff", yvar = "p_values"
    )
    protein_selected <- protein_tmp$name
  })
  protein_name_click <- reactive({
    protein_tmp <- nearPoints(volcano_df(), input$protein_click, maxpoints = 1)
    # protein_tmp<-brushedPoints(volcano_df(), input$protein_brush,
    # xvar = "diff", yvar = "p_values")
    protein_selected <- protein_tmp$name
  })


  ## Select rows dynamically

  makeReactiveBinding("brush")

  observeEvent(input$protein_brush, {
    output$contents <- DT::renderDataTable(
      {
        df <- data_result()[data_result()[["Gene Name"]] %in% protein_name_brush(), ]
        return(df)
      },
      options = list(scrollX = TRUE)
    )
  })

  observeEvent(input$resetPlot, {
    session$resetBrush("protein_brush")
    brush <<- NULL

    output$contents <- DT::renderDataTable(
      {
        df <- data_result()
        return(df)
      },
      options = list(
        scrollX = TRUE,
        autoWidth = TRUE,
        columnDefs = list(list(width = "400px", targets = c(-1)))
      )
    )
  })

  observeEvent(input$protein_click, {
    output$contents <- DT::renderDataTable(
      {
        df <- data_result()[data_result()[["Gene Name"]] %in% protein_name_click(), ]
        return(df)
      },
      options = list(
        scrollX = TRUE,
        autoWidth = TRUE,
        columnDefs = list(list(width = "400px", targets = c(-1)))
      )
    )
  })

  ## Render Result Plots
  output$pca_plot <- renderPlot({
    pca_input()
  })

  output$heatmap <- renderPlot({
    withProgress(
      message = "Heatmap rendering is in progress",
      detail = "Please wait for a while",
      value = 0,
      {
        for (i in 1:15) {
          incProgress(1 / 15)
          Sys.sleep(0.25)
        }
      }
    )
    heatmap_input()
  })

  output$volcano <- renderPlot({
    withProgress(
      message = "Volcano Plot calculations are in progress",
      detail = "Please wait for a while",
      value = 0,
      {
        for (i in 1:15) {
          incProgress(1 / 15)
          Sys.sleep(0.25)
        }
      }
    )
    if (is.null(input$contents_rows_selected) & is.null(input$protein_brush)) {
      volcano_input()
    } else if (!is.null(input$volcano_cntrst)) {
      volcano_input_selected()
    } # else close
  })

  output$protein_plot <- renderPlot({
    if (!is.null(input$contents_rows_selected)) {
      protein_input()
    }
  })


  ### QC Outputs
  output$sample_corr <- renderPlot({
    correlation_input()
  })

  output$sample_cvs <- renderPlot({
    cvs_input()
  })

  output$norm <- renderPlot({
    norm_input()
  })

  output$missval <- renderPlot({
    missval_input()
  })

  output$detect <- renderPlot({
    detect_input()
  })

  output$imputation <- renderPlot({
    imputation_input()
  })

  output$p_hist <- renderPlot({
    p_hist_input()
  })

  output$numbers <- renderPlot({
    numbers_input()
  })

  output$coverage <- renderPlot({
    coverage_input()
  })

  ## Enrichment Outputs
  output$spinner_go <- renderUI({
    req(input$go_analysis)
    shinycssloaders::withSpinner(plotOutput("go_enrichment"), color = "#3c8dbc")
  })

  output$go_enrichment <- renderPlot({
    Sys.sleep(2)
    go_input()$plot_go
  })

  output$spinner_pa <- renderUI({
    req(input$pathway_analysis)
    shinycssloaders::withSpinner(plotOutput("pathway_enrichment"), color = "#3c8dbc")
  })

  output$pathway_enrichment <- renderPlot({
    Sys.sleep(2)
    pathway_input()$plot_pa
  })

  ##### Download Functions
  datasetInput <- reactive({
    switch(input$dataset,
      "Results" = get_results_proteins(dep()),
      "Original_matrix" = unimputed_table(),
      # "significant_proteins" = get_results(dep()) %>%
      #   filter(significant) %>%
      #   select(-significant),
      "Imputed_matrix" = imputed_table(),
      "Full_dataset" = DEP::get_df_wide(dep())
    )
  })

  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    }, ## use = instead of <-
    content = function(file) {
      write.table(datasetInput(),
        file,
        col.names = TRUE,
        row.names = FALSE,
        sep = ","
      )
    }
  )

  ### === Cluster Download ==== ####

  individual_cluster <- reactive({
    cluster_number <- input$cluster_number
    cluster_all <- heatmap_cluster()[[2]]
    single_cluster <- cluster_all[names(cluster_all) == cluster_number] %>% unlist()
    df <- data_result()[single_cluster, ]
    return(df)
  })

  # output$text1 <- renderPrint({
  #   paste(individual_cluster())
  # })

  output$downloadCluster <- downloadHandler(
    filename = function() {
      paste("Cluster_info_", input$cluster_number, ".csv", sep = "")
    }, ## use = instead of <-
    content = function(file) {
      write.table(individual_cluster(),
        file,
        col.names = TRUE,
        row.names = FALSE,
        sep = ","
      )
    }
  )

  output$downloadVolcano <- downloadHandler(
    filename = function() {
      paste0("Volcano_", input$volcano_cntrst, ".pdf")
    },
    content = function(file) {
      pdf(file)
      if (is.null(input$contents_rows_selected) & is.null(input$protein_brush)) {
        print(volcano_input())
        dev.off()
      } else {
        # observeEvent(input$protein_brush,{
        #   print(p)
        # })
        print(volcano_input_selected())
        dev.off()
      }
    }
  )


  ## Protein plot download
  output$downloadProtein <- downloadHandler(
    filename = function() {
      paste0(input$type, ".pdf")
    },
    content = function(file) {
      pdf(file)
      print(protein_input())
      dev.off()
    }
  )

  ###### ==== DOWNLOAD GO TABLE ==== ####
  output$downloadGO <- downloadHandler(
    filename = function() {
      paste("GO_enrichment_", input$go_database, ".csv", sep = "")
    }, ## use = instead of <-
    content = function(file) {
      write.table(go_input()$go_result,
        file,
        col.names = TRUE,
        row.names = FALSE,
        sep = ","
      )
    }
  )

  ###### ==== DOWNLOAD PATHWAY TABLE ==== ####
  output$downloadPA <- downloadHandler(
    filename = function() {
      paste("Pathway_enrichment_", input$pathway_database, ".csv", sep = "")
    },
    ## use = instead of <-
    content = function(file) {
      write.table(pathway_input()$pa_result,
        file,
        col.names = TRUE,
        row.names = FALSE,
        sep = ","
      )
    }
  )

  output$download_hm_svg <- downloadHandler(
    filename = function() {
      "heatmap.svg"
    },
    ## use = instead of <-
    content = function(file) {
      # heatmap_plot<-DEP::plot_heatmap(dep(),"centered", k=6, indicate = "condition")
      svg(file)
      print(heatmap_input())
      dev.off()
    }
  )

  ##### ===== Download Report =====#####
  output$downloadReport <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "LFQ-Analyst_report.pdf",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "LFQ_report.Rmd")
      file.copy("LFQ_report.Rmd", tempReport, overwrite = TRUE)

      sig_proteins <- dep() %>%
        .[SummarizedExperiment::rowData(.)$significant, ] %>%
        nrow()

      tested_contrasts <- gsub(
        "_p.adj", "",
        colnames(SummarizedExperiment::rowData(dep()))[grep(
          "p.adj",
          colnames(SummarizedExperiment::rowData(dep()))
        )]
      )
      pg_width <- ncol(imputed_data()) / 2.5
      # Set up parameters to pass to Rmd document
      params <- list(
        data = processed_data,
        alpha = input$p,
        lfc = input$lfc,
        num_signif = sig_proteins,
        pg_width = pg_width,
        tested_contrasts = tested_contrasts,
        numbers_input = numbers_input,
        detect_input = detect_input,
        imputation_input = imputation_input,
        missval_input = missval_input,
        p_hist_input = p_hist_input,
        pca_input = pca_input,
        coverage_input = coverage_input,
        correlation_input = correlation_input,
        heatmap_input = heatmap_input,
        cvs_input = cvs_input,
        dep = dep
      )

      # Knit the document, passing in the `params` list
      rmarkdown::render(tempReport,
        output_file = file,
        params = params,
        envir = new.env(parent = globalenv())
      )
    }
  )

  ###### ==== DOWNLOAD QC plots svg ==== ####

  output$download_pca_svg <- downloadHandler(
    filename = function() {
      "PCA_plot.svg"
    },
    content = function(file) {
      svg(file)
      print(pca_input())
      dev.off()
    }
  )

  output$download_corr_svg <- downloadHandler(
    filename = function() {
      "Correlation_plot.svg"
    },
    content = function(file) {
      svg(file)
      print(correlation_input())
      dev.off()
    }
  )

  output$download_cvs_svg <- downloadHandler(
    filename = function() {
      "Sample_CV.svg"
    },
    content = function(file) {
      svg(file)
      print(cvs_input())
      dev.off()
    }
  )

  output$download_num_svg <- downloadHandler(
    filename = function() {
      "Proteins_plot.svg"
    },
    content = function(file) {
      svg(file)
      print(numbers_input())
      dev.off()
    }
  )

  output$download_cov_svg <- downloadHandler(
    filename = function() {
      "Coverage_plot.svg"
    },
    content = function(file) {
      svg(file)
      print(coverage_input())
      dev.off()
    }
  )

  output$download_norm_svg <- downloadHandler(
    filename = function() {
      "Normalization_plot.svg"
    },
    content = function(file) {
      svg(file)
      print(norm_input())
      dev.off()
    }
  )

  output$download_missval_svg <- downloadHandler(
    filename = function() {
      "Missing_value_heatmap.svg"
    },
    content = function(file) {
      svg(file)
      print(missval_input())
      dev.off()
    }
  )

  output$download_imp_svg <- downloadHandler(
    filename = function() {
      "Imputation_plot.svg"
    },
    content = function(file) {
      svg(file)
      print(imputation_input())
      dev.off()
    }
  )
}
