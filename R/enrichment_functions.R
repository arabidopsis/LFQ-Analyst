# Term <- ""
# Overlap <- ""
# IN <- ""
# n <- ""
# bg_OUT <- ""
# bg_IN <- ""
# name <- ""

`%>%` <- magrittr::`%>%`

enrichr_mod <- function(genes, databases = NULL) {
  httr::set_config(httr::config(ssl_verifypeer = 0L))
  cat("Uploading data to Enrichr... ")
  if (is.vector(genes) && !all(genes == "") && length(genes) != 0) {
    httr::POST(
      url = "http://maayanlab.cloud/Enrichr/enrich",
      body = list(list = paste(genes, collapse = "\n"))
    )
  } else if (is.data.frame(genes)) {
    httr::POST(
      url = "http://maayanlab.cloud/Enrichr/enrich",
      body = list(list = paste(paste(genes[, 1], genes[, 2], sep = ","),
        collapse = "\n"
      ))
    )
  } else {
    warning("genes must be a non-empty vector of gene names or a dataframe with genes and score.")
  }
  httr::GET(url = "http://maayanlab.cloud/Enrichr/share")
  cat("Done.\n")
  dbs <- as.list(databases)
  dfSAF <- options()$stringsAsFactors
  options(stringsAsFactors = FALSE)
  result <- lapply(dbs, function(x) {
    cat("  Querying ", x, "... ", sep = "")
    r <- httr::GET(
      url = "http://maayanlab.cloud/Enrichr/export",
      query = list(file = "API", backgroundType = x)
    )
    r <- gsub("&#39;", "'", intToUtf8(r$content))
    tc <- textConnection(r)
    r <- read.table(tc, sep = "\t", header = TRUE, quote = "", comment.char = "")
    close(tc)
    cat("Done.\n")
    r
  })
  options(stringsAsFactors = dfSAF)
  cat("Parsing results... ")
  names(result) <- dbs
  cat("Done.\n")
  return(result)
}


###### ========= Test_gsea new


test_gsea_mod <- function(dep,
                          databases,
                          contrasts = TRUE) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(
    inherits(dep, "SummarizedExperiment"),
    is.character(databases),
    is.logical(contrasts),
    length(contrasts) == 1
  )


  row_data <- SummarizedExperiment::rowData(dep, use.names = FALSE)
  # Show error if inputs do not contain required columns
  if (any(!c("name", "ID") %in% colnames(row_data))) {
    stop("'name' and/or 'ID' columns are not present in '",
      deparse(substitute(dep)),
      "'\nRun make_unique() and make_se() to obtain the required columns",
      call. = FALSE
    )
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '",
      deparse(substitute(dep)),
      "'\nRun test_diff() to obtain the required columns",
      call. = FALSE
    )
  }



  # Run background list
  message("Background")
  background <- gsub("[.].*", "", row_data$name)
  background_enriched <- enrichr_mod(background, databases)
  df_background <- NULL
  for (database in databases) {
    temp <- background_enriched[database][[1]] %>%
      dplyr::mutate(var = database)
    df_background <- rbind(df_background, temp)
  }
  df_background$contrast <- "background"
  df_background$n <- length(background)

  OUT <- df_background %>%
    dplyr::mutate(
      bg_IN = as.numeric(gsub("/.*", "", Overlap)),
      bg_OUT = n - bg_IN
    ) %>%
    dplyr::select(Term, bg_IN, bg_OUT)

  if (contrasts) {
    # Get gene symbols
    df <- row_data %>%
      as.data.frame() %>%
      dplyr::select(name, dplyr::ends_with("_significant")) %>%
      dplyr::mutate(name = gsub("[.].*", "", name))

    # Run enrichR for every contrast
    df_enrich <- NULL
    for (contrast in colnames(df[2:ncol(df)])) {
      message(gsub("_significant", "", contrast))
      significant <- df[df[[contrast]], ]
      genes <- significant$name
      enriched <- enrichr_mod(genes, databases)

      # Tidy output
      contrast_enrich <- NULL
      for (database in databases) {
        temp <- enriched[database][[1]] %>%
          dplyr::mutate(var = database)
        contrast_enrich <- rbind(contrast_enrich, temp)
      }

      if (nrow(contrast_enrich) != 0) {
        contrast_enrich$contrast <- contrast
        contrast_enrich$n <- length(genes)

        # Background correction
        cat("Background correction... ")
        contrast_enrich <- contrast_enrich %>%
          dplyr::mutate(
            IN = as.numeric(gsub("/.*", "", Overlap)),
            OUT = n - IN
          ) %>%
          dplyr::select(-n) %>%
          dplyr::left_join(OUT, by = "Term") %>%
          dplyr::mutate(log_odds = log2((IN * bg_OUT) / (OUT * bg_IN)))
        cat("Done.")

        df_enrich <- rbind(df_enrich, contrast_enrich) %>%
          dplyr::mutate(contrast = gsub("_significant", "", contrast))
      } else {
        cat("No enough significant genes for enrichment analysis")
      }
    }
  } else {
    # Get gene symbols
    significant <- row_data %>%
      as.data.frame() %>%
      dplyr::select(name, significant) %>%
      filter(significant) %>%
      dplyr::mutate(name = gsub("[.].*", "", name))

    # Run enrichR
    genes <- significant$name
    enriched <- enrichr_mod(genes, databases)

    # Tidy output
    df_enrich <- NULL
    for (database in databases) {
      temp <- enriched[database][[1]] %>%
        dplyr::mutate(var = database)
      df_enrich <- rbind(df_enrich, temp)
    }
    df_enrich$contrast <- "significant"
    df_enrich$n <- length(genes)

    # Background correction
    cat("Background correction... ")
    df_enrich <- df_enrich %>%
      dplyr::mutate(
        IN = as.numeric(gsub("/.*", "", Overlap)),
        OUT = n - IN
      ) %>%
      dplyr::select(-n) %>%
      dplyr::left_join(OUT, by = "Term") %>%
      dplyr::mutate(log_odds = log2((IN * bg_OUT) / (OUT * bg_IN)))
    cat("Done.")
  }

  return(df_enrich)
}


# gene_names_true<-read_table("R/gene_names.txt",col_names = F)
