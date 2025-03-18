# Define server logic to read selected file ----
library("shiny", quietly = TRUE)


name_space <- shiny::NS("lfq")

significantBox <- function(row_data, max_frac = 0.2) {
  num_signif <- nrow(row_data[row_data$significant, ])
  num_total <- nrow(row_data)

  frac <- num_signif / num_total

  value <- paste0(signif(frac * 100, digits = 3), "%")
  title <- paste0(
    num_signif,
    " out of ",
    num_total,
    " proteins differentially expressed"
  )
  if (frac > max_frac) {
    title <- paste(title, "(too large!)")
    icon <- "minus"
    theme <- "text-warning"
  } else if (frac == 0) {
    icon <- "thumbs-down"
    theme <- "text-danger"
  } else {
    icon <- "thumbs-up"
    theme <- "text-success"
  }
  bslib::value_box(
    title = title,
    value = value,
    showcase = shiny::icon(icon),
    theme = theme,
    height = "100%"
  )
}

server_bg <- function(input, output, session) {
  #  Show elements on clicking Start analysis button
  observeEvent(input$analyze,
    {
      if (input$analyze == 0) {
        return()
      }
      shinyjs::hide("quickstart_info")
      shinyjs::show("analysis_id")

      shinyalert::shinyalert("In Progress!", "Data analysis has started, wait until table and plots
                appear on the screen",
        type = "info",
        closeOnClickOutside = TRUE,
        closeOnEsc = TRUE,
        timer = 5000
      ) # timer in miliseconds (10 sec)
    },
    ignoreNULL = TRUE
  )

  maxquant_data_input <- eventReactive(input$analyze,
    {
      inFile <- input$maxquant_file
      if (is.null(inFile)) {
        return(NULL)
      }
      temp_data <- read.table(inFile$datapath,
        header = TRUE,
        fill = TRUE, # to fill any missing data
        sep = "\t",
        quote = ""
      )
      maxquant_input_test(temp_data)
      return(temp_data)
    },
    ignoreNULL = TRUE
  )

  exp_design_input <- eventReactive(input$analyze,
    {
      inFile <- input$exp_design_file
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
      temp_df$condition <- gsub("[^[:alnum:]|_]+", "_", temp_df$condition) # auto fix special characters
      return(temp_df)
    },
    ignoreNULL = TRUE
  )

  # these need to be global because of LFQ_report.Rmd

  processed_data <- reactive({
    if (any(grepl("+", maxquant_data_input()$Reverse))) {
      filtered_data <- dplyr::filter(maxquant_data_input(), Reverse != "+")
    } else {
      filtered_data <- maxquant_data_input()
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
    test_match_lfq_column_design(data_unique, lfq_columns, exp_design_input())
    data_se <- DEP::make_se(data_unique, lfq_columns, exp_design_input())


    exp_df <- dplyr::count(exp_design_input(), condition)
    exp_df <- dplyr::mutate(exp_df, thr = lapply(exp_df$n, threshold_detect)) # function:threshold_detect
    condition_list <- exp_df$condition

    data_filtered <- filter_missval_new(data_se, condition_list, exp_df)
    return(data_filtered)
  })

  normalised_data <- reactive({
    DEP::normalize_vsn(processed_data())
  })

  imputed_data <- reactive({
    DEP::impute(processed_data(), input$imputation)
  })

  diff_all <- reactive({
    if (input$fdr_correction == "BH") {
      diff <- test_limma(imputed_data(), type = "all", paired = input$paired)
    } else {
      diff <- DEP::test_diff(imputed_data(), type = "all")
    }
    diff
  })

  dep <- reactive({
    DEP::add_rejections(diff_all(), alpha = input$p_value, lfc = input$log_fold_change)
  })

  # comparisons <- reactive({
  #   temp <- capture.output(DEP::test_diff(imputed_data(), type = "all"), type = "message")
  #   gsub(".*: ", "", temp)
  #   ## Split conditions into character vector
  #   unlist(strsplit(temp, ","))
  #   ## Remove leading and trailing spaces
  #   trimws(temp)
  # })

  comparisons2 <- reactive({
    df <- SummarizedExperiment::rowData(dep())
    cols <- grep("_significant$", colnames(df))
    gsub("_significant", "", colnames(df)[cols])
  })

  ### Heatmap Differentially expressed proteins
  heatmap_cluster <- reactive({
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
      if (!input$p_adj) {
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
        label_size = input$fontsize,
        add_names = input$check_names,
        adjusted = input$p_adj
      )

      p + ggplot2::geom_point(data = df_protein, ggplot2::aes(x, y), color = "maroon", size = 3) +
        ggrepel::geom_text_repel(
          data = df_protein,
          ggplot2::aes(x, y, label = name),
          size = 4,
          box.padding = grid::unit(0.1, "lines"),
          point.padding = grid::unit(0.1, "lines"),
          segment.size = 0.5
        ) ## use the dataframe to plot points
    }
  })

  ## QC Inputs
  pca_input <- reactive({
    num_total <- nrow(dep())
    if (num_total <= 500) {
      if (length(levels(as.factor(SummarizedExperiment::colData(dep())$replicate))) <= 6) {
        pca_plot <- DEP::plot_pca(dep(), n = num_total(), point_size = 4)
        pca_plot <- pca_plot + ggplot2::labs(title = "PCA Plot")
        return(pca_plot)
      } else {
        pca_plot <- DEP::plot_pca(dep(), n = num_total(), point_size = 4, indicate = "condition")
        pca_plot <- pca_plot + ggplot2::labs(title = "PCA Plot")
        return(pca_plot)
      }
    } else {
      if (length(levels(as.factor(SummarizedExperiment::colData(dep())$replicate))) <= 6) {
        pca_plot <- DEP::plot_pca(dep(), point_size = 4)
        pca_plot <- pca_plot + ggplot2::labs(title = "PCA Plot")
        return(pca_plot)
      } else {
        # pca_label<-SummarizedExperiment::colData(dep())$replicate
        pca_plot <- DEP::plot_pca(dep(), point_size = 4, indicate = "condition")
        # pca_plot<-pca_plot + geom_point()
        pca_plot <- pca_plot + ggrepel::geom_text_repel(ggplot2::aes(label = factor(rowname)),
          size = 4,
          box.padding = grid::unit(0.1, "lines"),
          point.padding = grid::unit(0.1, "lines"),
          segment.size = 0.5
        )
        pca_plot <- pca_plot + ggplot2::labs(title = "PCA Plot")
        return(pca_plot)
      }
    }
  })

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
    nrow(dep())
  })

  data_result <- reactive({
    get_results_proteins(dep())
    # get_results(dep())
  })

  protein_name_brush <- reactive({
    # protein_tmp<-nearPoints(volcano_df(), input$protein_click, maxpoints = 1)
    protein_tmp <- brushedPoints(volcano_df(), input$protein_brush,
      xvar = "diff", yvar = "p_values"
    )
    protein_tmp$name
  })

  protein_name_click <- reactive({
    protein_tmp <- nearPoints(volcano_df(), input$protein_click, maxpoints = 1)
    # protein_tmp<-brushedPoints(volcano_df(), input$protein_brush,
    # xvar = "diff", yvar = "p_values")
    protein_tmp$name
  })

  #### ==== table panel ==== ####
  table_panel <- function() {
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
    ## Select rows dynamically

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
      # brush <<- NULL

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



    output$heatmap_plot <- renderPlot({
      # withProgress(
      #   message = "Heatmap rendering is in progress",
      #   detail = "Please wait for a while",
      #   value = 0,
      #   session = session,
      #   {
      #     for (i in 1:15) {
      #       incProgress(1 / 15)
      #       Sys.sleep(0.25)
      #     }
      #   }
      # )
      heatmap_input()
    })
  }

  table_panel()

  #### ==== top row panel ==== ####
  top_row_panel <- function() {

    unimputed_table <- reactive({
      temp <- SummarizedExperiment::assay(processed_data())
      temp1 <- 2^(temp)
      colnames(temp1) <- paste(colnames(temp1), "original_intensity", sep = "_")
      temp1 <- cbind(ProteinID = rownames(temp1), temp1)
      # temp1$ProteinID<-rownames(temp1)
      return(as.data.frame(temp1))
    })

    imputed_table <- reactive({
      temp <- SummarizedExperiment::assay(imputed_data())
      # tibble::rownames_to_column(temp,var = "ProteinID")
      temp1 <- 2^(temp)
      colnames(temp1) <- paste(colnames(temp1), "imputed_intensity", sep = "_")
      temp1 <- cbind(ProteinID = rownames(temp1), temp1) # temp1$ProteinID<-rownames(temp1)
      return(as.data.frame(temp1))
    })

    output$significantBox <- renderUI({
      significantBox(SummarizedExperiment::rowData(dep()))
    })

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
          alpha = input$p_value,
          log_fold_change = input$log_fold_change,
          num_signif = sig_proteins,
          pg_width = pg_width,
          tested_contrasts = tested_contrasts,
          numbers_input = numbers_input,
          detect_input = detect_input,
          imputation_input = imputation_input,
          missval_input = missval_input,
          pca_input = pca_input,
          coverage_input = coverage_input,
          correlation_input = correlation_input,
          heatmap_input = heatmap_input,
          cvs_input = cvs_input,
          dep = dep
        )

        # Knit the document, passing in the `params` list
        tryCatch(
          {
            rmarkdown::render(tempReport,
              output_file = file,
              params = params,
              envir = new.env(parent = globalenv())
            )
          },
          finally = {
            file.remove(tempReport)
          }
        )
      }
    )
  }

  top_row_panel()

  #### ==== volcano panel ==== ####
  volcano_panel <- function() {
    volcano_input <- reactive({
      if (!is.null(input$volcano_cntrst)) {
        plot_volcano_new(
          dep(),
          input$volcano_cntrst,
          label_size = input$fontsize,
          add_names = input$check_names,
          adjusted = input$p_adj
        )
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
        plot_protein(dep(), protein_selected, input$protein_plot_type)
      } else {
        pp <- plot_protein(dep(), protein_selected, input$protein_plot_type)
        pp + ggplot2::scale_color_brewer(palette = "Paired")
      }
    })

    individual_cluster <- reactive({
      cluster_number <- input$cluster_number
      cluster_all <- heatmap_cluster()[[2]]
      single_cluster <- cluster_all[names(cluster_all) == cluster_number] %>% unlist()
      df <- data_result()[single_cluster, ]
      return(df)
    })

    #### ======= Render Functions
    output$volcano_cntrst_placeholder <- renderUI({
      selectizeInput(name_space("volcano_cntrst"),
        "Comparison",
        choices = comparisons2(),
        options = list(dropdownParent = "body")
      )
    })

    output$protein_plot <- renderPlot({
      if (!is.null(input$contents_rows_selected)) {
        protein_input()
      }
    })

    output$download_hm_svg <- downloadHandler(
      filename = function() {
        "heatmap.svg"
      },
      ## use = instead of <-
      content = function(file) {
        svg(file)
        print(heatmap_input())
        dev.off()
      }
    )
    output$volcano_plot <- renderPlot({
      # withProgress(
      #   message = "Volcano Plot calculations are in progress",
      #   detail = "Please wait for a while",
      #   value = 0,
      #   session = session,
      #   {
      #     for (i in 1:15) {
      #       incProgress(1 / 15)
      #       Sys.sleep(0.25)
      #     }
      #   }
      # )
      if (is.null(input$contents_rows_selected) && is.null(input$protein_brush)) {
        volcano_input()
      } else if (!is.null(input$volcano_cntrst)) {
        volcano_input_selected()
      } # else close
    })

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
        paste0(input$protein_plot_type, ".pdf")
      },
      content = function(file) {
        pdf(file)
        print(protein_input())
        dev.off()
      }
    )
  }

  volcano_panel()

  #### ==== enrichment panel ==== ####
  enrichment_panel <- function() {
    go_input <- eventReactive(input$go_analysis, {
      withProgress(
        message = "Gene ontology enrichment is in progress",
        detail = "Please wait for a while",
        value = 0,
        session = session,
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
        ) + ggplot2::aes(stringr::str_wrap(Term, 60)) +
          ggplot2::xlab(NULL)
        go_list <- list("go_result" = go_results, "plot_go" = plot_go)
        return(go_list)
      }
    })
    pathway_input <- eventReactive(input$pathway_analysis, {
      shiny::withProgress(
        message = "Pathway Analysis is running....",
        detail = "Please wait for a while",
        value = 0,
        {
          for (i in 1:15) {
            shiny::incProgress(1 / 15)
            Sys.sleep(0.25)
          }
        }
      )
      enrichment_output_test(dep(), input$pathway_database)
      pathway_results <- test_gsea_mod(dep(), databases = input$pathway_database, contrasts = TRUE)
      null_enrichment_test(pathway_results)
      plot_pathway <- plot_enrichment(pathway_results,
        number = 5, alpha = 0.05, contrasts = input$contrast_1,
        databases = input$pathway_database, nrow = 3, term_size = 8
      ) + ggplot2::aes(stringr::str_wrap(Term, 30)) +
        ggplot2::xlab(NULL)
      pathway_list <- list("pa_result" = pathway_results, "plot_pa" = plot_pathway)
      return(pathway_list)
    })
    ## Enrichment Outputs


    output$contrast_placeholder <- renderUI({
      selectizeInput(name_space("contrast"),
        "Comparison",
        choices = comparisons2(),
        options = list(dropdownParent = "body")
      )
    })

    output$contrast_1_placeholder <- renderUI({
      selectizeInput(name_space("contrast_1"),
        "Comparison",
        choices = comparisons2(),
        options = list(dropdownParent = "body")
      )
    })
    output$spinner_go <- renderUI({
      req(input$go_analysis)
      shinycssloaders::withSpinner(plotOutput(name_space("go_enrichment")), color = "#3c8dbc")
    })

    output$go_enrichment <- renderPlot({
      Sys.sleep(2)
      go_input()$plot_go
    })

    output$spinner_pa <- renderUI({
      req(input$pathway_analysis)
      shinycssloaders::withSpinner(plotOutput(name_space("pathway_enrichment")), color = "#3c8dbc")
    })

    output$pathway_enrichment <- renderPlot({
      Sys.sleep(2)
      pathway_input()$plot_pa
    })

    output$spinner_go <- renderUI({
      req(input$go_analysis)
      shinycssloaders::withSpinner(plotOutput(name_space("go_enrichment")), color = "#3c8dbc")
    })

    output$go_enrichment <- renderPlot({
      Sys.sleep(2)
      go_input()$plot_go
    })

    output$spinner_pa <- renderUI({
      req(input$pathway_analysis)
      shinycssloaders::withSpinner(plotOutput(name_space("pathway_enrichment")), color = "#3c8dbc")
    })

    output$pathway_enrichment <- renderPlot({
      Sys.sleep(2)
      pathway_input()$plot_pa
    })

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
  }

  enrichment_panel()

  #### ==== QC panel ==== ####
  qc_panel <- function() {
    ### QC Outputs
    output$pca_plot <- renderPlot({
      pca_input()
    })

    output$sample_corr_plot <- renderPlot({
      correlation_input()
    })

    output$sample_cvs_plot <- renderPlot({
      cvs_input()
    })

    output$normalization_plot <- renderPlot({
      norm_input()
    })

    output$missval_plot <- renderPlot({
      missval_input()
    })

    output$detect <- renderPlot({
      detect_input()
    })

    output$imputation_plot <- renderPlot({
      imputation_input()
    })

    output$numbers_plot <- renderPlot({
      numbers_input()
    })

    output$coverage_plot <- renderPlot({
      coverage_input()
    })

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

  qc_panel()
}

server <- function(input, output, session) {
  options(shiny.maxRequestSize = 100 * 1024^2) ## Set maximum upload size to 100MB
  moduleServer("lfq", server_bg)
  # server_bg(input, output, session)
}
