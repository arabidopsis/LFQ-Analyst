# library(ggplot2)

# ID <- ""
# condition <- ""
# cvs <- ""
# Intensity <- ""
# condition_median <- ""
# cluster <- ""
# Reverse <- ""
# Only.identified.by.site <- ""
# rowname <- ""
# index <- ""
# logFC <- ""
# CI.L <- ""
# CI.R <- ""
# P.Value <- ""
# adj.P.Val <- ""
# comparison <- ""
# variable <- ""
# val <- ""
# name <- ""
# significant <- ""
# contrast <- ""
# Adjusted.P.value <- ""
# Term <- ""
# miss_val <- ""
# Protein.names <- ""
# value <- ""
# log_odds <- ""


`%>%` <- magrittr::`%>%`


coef_variation <- function(x) {
  sd(x) / mean(x)
}

#### Plot CVs

#' @export
plot_cvs <- function(se) {
  ## backtransform data
  untransformed_intensity <- 2^(SummarizedExperiment::assay(se))
  exp_design <- SummarizedExperiment::colData(se)

  ### merge untransformed to exp design and calculate cvs

  cvs_group <- untransformed_intensity %>%
    data.frame() %>%
    tibble::rownames_to_column() %>%
    tidyr::gather("ID", "Intensity", -rowname) %>%
    dplyr::left_join(data.frame(exp_design), by = "ID") %>%
    dplyr::group_by(rowname, condition) %>%
    dplyr::summarise(cvs = coef_variation(Intensity)) %>%
    dplyr::group_by(condition) %>%
    dplyr::mutate(condition_median = median(cvs))

  p1 <- ggplot2::ggplot(cvs_group, ggplot2::aes(cvs, color = condition, fill = condition)) +
    ggplot2::geom_histogram(alpha = .5, bins = 20, show.legend = FALSE) +
    ggplot2::facet_wrap(~condition) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = condition_median, group = condition),
      color = "grey40",
      linetype = "dashed"
    ) +
    ggplot2::labs(title = "Sample Coefficient of Variation", x = "Coefficient of Variation", y = "Count") +
    DEP::theme_DEP2() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  p1 + ggplot2::geom_text(
    ggplot2::aes(
      x = max(cvs_group$cvs) - 0.6,
      y = max(ggplot2::ggplot_build(p1)$data[[1]]$ymax * 1.1),
      label = paste0("Median =", round(condition_median, 2) * 100, "%", by = "")
    ),
    show.legend = FALSE, size = 4
  )
}


#### Get individual clusters from heatmap
#' @export
get_cluster_heatmap <- function(dep, type = c("contrast", "centered"),
                                kmeans = FALSE, k = 6,
                                col_limit = 6, indicate = NULL,
                                clustering_distance = c(
                                  "euclidean", "maximum", "manhattan", "canberra",
                                  "binary", "minkowski", "pearson", "spearman", "kendall", "gower"
                                ),
                                row_font_size = 6, col_font_size = 10, plot = TRUE, ...) {
  # Show error if inputs are not the required classes
  if (is.integer(k)) k <- as.numeric(k)
  if (is.integer(col_limit)) col_limit <- as.numeric(col_limit)
  if (is.integer(row_font_size)) row_font_size <- as.numeric(row_font_size)
  if (is.integer(col_font_size)) col_font_size <- as.numeric(col_font_size)
  assertthat::assert_that(
    inherits(dep, "SummarizedExperiment"),
    is.character(type),
    is.logical(kmeans),
    is.numeric(k),
    length(k) == 1,
    is.numeric(col_limit),
    length(col_limit) == 1,
    is.numeric(row_font_size),
    length(row_font_size) == 1,
    is.numeric(col_font_size),
    length(col_font_size) == 1,
    is.logical(plot),
    length(plot) == 1
  )

  # Show error if inputs do not contain required columns
  type <- match.arg(type)
  clustering_distance <- match.arg(clustering_distance)

  # Extract row and col data
  row_data <- SummarizedExperiment::rowData(dep)
  col_data <- SummarizedExperiment::colData(dep) %>%
    as.data.frame()

  # Show error if inputs do not contain required columns
  if (any(!c("label", "condition", "replicate") %in% colnames(col_data))) {
    stop(
      paste0(
        "'label', 'condition' and/or 'replicate' columns are not present in '",
        deparse(substitute(dep)), "'"
      ),
      call. = FALSE
    )
  }
  if (length(grep("_diff", colnames(row_data))) < 1) {
    stop(
      paste0(
        "'[contrast]_diff' columns are not present in '",
        deparse(substitute(dep)),
        "'.\nRun test_diff() to obtain the required columns."
      ),
      call. = FALSE
    )
  }
  if (!"significant" %in% colnames(row_data)) {
    stop(
      paste0(
        "'significant' column is not present in '",
        deparse(substitute(dep)),
        "'.\nRun add_rejections() to obtain the required column."
      ),
      call. = FALSE
    )
  }

  # Heatmap annotation
  if (!is.null(indicate) && type == "contrast") {
    warning("Heatmap annotation only applicable for type = 'centered'",
      call. = FALSE
    )
  }
  if (!is.null(indicate) && type == "centered") {
    ha1 <- get_annotation(dep, indicate)
  } else {
    ha1 <- NULL
  }

  # Filter for significant proteins only
  filtered <- dep[row_data$significant, ]

  # Check for missing values
  if (any(is.na(SummarizedExperiment::assay(filtered)))) {
    warning("Missing values in '", deparse(substitute(dep)), "'. ",
      "Using clustering_distance = 'gower'",
      call. = FALSE
    )
    clustering_distance <- "gower"
    obs_NA <- TRUE
  } else {
    obs_NA <- FALSE
  }

  # Get centered intensity values ('centered')
  if (type == "centered") {
    SummarizedExperiment::rowData(filtered)$mean <- rowMeans(SummarizedExperiment::assay(filtered), na.rm = TRUE)
    df <- SummarizedExperiment::assay(filtered) - SummarizedExperiment::rowData(filtered)$mean
  }
  # Get contrast fold changes ('contrast')
  if (type == "contrast") {
    df <- SummarizedExperiment::rowData(filtered) %>%
      data.frame() %>%
      tibble::column_to_rownames(var = "name") %>%
      dplyr::select(tidyselect::ends_with("_diff"))
    colnames(df) <-
      gsub("_diff", "", colnames(df)) %>%
      gsub("_vs_", " vs ", .)
  }

  # Facultative kmeans clustering
  if (kmeans && obs_NA) {
    warning("Cannot perform kmeans clustering with missing values",
      call. = FALSE
    )
    kmeans <- FALSE
  }
  if (kmeans && !obs_NA) {
    set.seed(1)
    df_kmeans <- kmeans(df, k)
    if (type == "centered") {
      # Order the k-means clusters according to the maximum fold change
      # in all samples averaged over the proteins in the cluster
      order <- data.frame(df) %>%
        cbind(cluster = df_kmeans$cluster) %>%
        dplyr::mutate(row = apply(.[, seq_len(ncol(.) - 1)], 1, function(x) max(x))) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarize(index = sum(row) / dplyr::n()) %>%
        dplyr::arrange(dplyr::desc(index)) %>%
        dplyr::pull(cluster) %>%
        match(seq_len(k), .)
      df_kmeans$cluster <- order[df_kmeans$cluster]
    }
    if (type == "contrast") {
      # Order the k-means clusters according to their average fold change
      order <- cbind(df, cluster = df_kmeans$cluster) %>%
        tidyr::gather(condition, diff, -cluster) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarize(row = mean(diff)) %>%
        dplyr::arrange(dplyr::desc(row)) %>%
        dplyr::pull(cluster) %>%
        match(seq_len(k), .)
      df_kmeans$cluster <- order[df_kmeans$cluster]
    }
  }

  if (ncol(df) == 1) {
    col_clust <- FALSE
  } else {
    col_clust <- TRUE
  }
  if (nrow(df) == 1) {
    row_clust <- FALSE
  } else {
    row_clust <- TRUE
  }
  if (clustering_distance == "gower") {
    clustering_distance <- function(x) {
      dist <- cluster::daisy(x, metric = "gower")
      dist[is.na(dist)] <- max(dist, na.rm = TRUE)
      dist
    }
  }

  # Legend info
  legend <- ifelse(type == "contrast",
    "log2 Fold change",
    "log2 Centered intensity"
  )

  # Heatmap
  ht1 <- ComplexHeatmap::Heatmap(df,
    col = circlize::colorRamp2(
      seq(-col_limit, col_limit, (col_limit / 5)),
      rev(RColorBrewer::brewer.pal(11, "RdBu"))
    ),
    split = if (kmeans) {
      df_kmeans$cluster
    } else {
      NULL
    },
    cluster_rows = col_clust,
    cluster_columns = row_clust,
    row_names_side = "left",
    column_names_side = "top",
    clustering_distance_rows = clustering_distance,
    clustering_distance_columns = clustering_distance,
    heatmap_legend_param = list(
      color_bar = "continuous",
      legend_direction = "horizontal",
      legend_width = grid::unit(5, "cm"),
      title_position = "lefttop"
    ),
    name = legend,
    row_names_gp = grid::gpar(fontsize = row_font_size),
    column_names_gp = grid::gpar(fontsize = col_font_size),
    top_annotation = ha1,
    ...
  )
  # return (row_order(ht1))
  # Return data.frame
  p <- ComplexHeatmap::draw(ht1, heatmap_legend_side = "top")
  row_clusters <- ComplexHeatmap::row_order(ht1)
  # mat<-as.matrix(df)

  # for (i in 1:length(row_clusters)){
  #   if (i==1){
  #     clu <-t(t(row.names(ht1[row_clusters[[i]],])))
  #     out <-cbind (clu, paste("cluster", i, sep=""))
  #     colnames(out)<- c("ProteinID", "Cluster")
  #   }
  #   else{
  #     clu <- t(t(row.names(ht1[row_clusters[[i]],])))
  #     clu <- cbind(clu, paste("cluster", i, sep = ""))
  #     out <- cbind(out, clu)
  #   }
  # }
  heatmap_list <- list(p, row_clusters)
  return(heatmap_list)
}

# Internal function to get ComplexHeatmap::HeatmapAnnotation object
get_annotation <- function(dep, indicate) {
  assertthat::assert_that(
    inherits(dep, "SummarizedExperiment"),
    is.character(indicate)
  )

  # Check indicate columns
  col_data <- SummarizedExperiment::colData(dep) %>%
    as.data.frame()
  columns <- colnames(col_data)
  if (all(!indicate %in% columns)) {
    stop("'",
      paste0(indicate, collapse = "' and/or '"),
      "' column(s) is/are not present in ",
      deparse(substitute(dep)),
      ".\nValid columns are: '",
      paste(columns, collapse = "', '"),
      "'.",
      call. = FALSE
    )
  }
  if (any(!indicate %in% columns)) {
    indicate <- indicate[indicate %in% columns]
    warning(
      "Only used the following indicate column(s): '",
      paste0(indicate, collapse = "', '"),
      "'"
    )
  }

  # Get annotation
  anno <- dplyr::select(col_data, indicate)

  # Annotation color
  names <- colnames(anno)
  anno_col <- vector(mode = "list", length = length(names))
  names(anno_col) <- names
  for (i in names) {
    var <- anno[[i]] %>%
      unique() %>%
      sort()
    if (length(var) == 1) {
      cols <- c("black")
    }
    if (length(var) == 2) {
      cols <- c("orangered", "cornflowerblue")
    }
    if (length(var) < 7 && length(var) > 2) {
      cols <- RColorBrewer::brewer.pal(length(var), "Pastel1")
    }
    if (length(var) > 7) {
      cols <- RColorBrewer::brewer.pal(length(var), "Set3")
    }
    names(cols) <- var
    anno_col[[i]] <- cols
  }

  # HeatmapAnnotation object
  ComplexHeatmap::HeatmapAnnotation(
    df = anno,
    col = anno_col,
    show_annotation_name = TRUE
  )
}


#### ===== limma BH FDR ===== #####

#' @export
test_limma <- function(se, type = c("control", "all", "manual"),
                       control = NULL, test = NULL,
                       design_formula = formula(~ 0 + condition),
                       paired = FALSE) {
  # require("dplyr", "tidyr", "purrr")

  # Show error if inputs are not the required classes
  assertthat::assert_that(
    inherits(se, "SummarizedExperiment"),
    is.character(type),
    class(design_formula) == "formula"
  )
  if (paired == FALSE) {
    design_formula <- design_formula
  } else {
    design_formula <- formula(~ 0 + condition + replicate)
  }


  # Show error if inputs do not contain required columns
  type <- match.arg(type)

  col_data <- SummarizedExperiment::colData(se)
  raw <- SummarizedExperiment::assay(se)

  if (any(!c("name", "ID") %in% colnames(SummarizedExperiment::rowData(se)))) {
    stop("'name' and/or 'ID' columns are not present in '",
      deparse(substitute(se)),
      "'\nRun make_unique() and make_se() to obtain the required columns",
      call. = FALSE
    )
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(col_data))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
      deparse(substitute(se)),
      "'\nRun make_se() or make_se_parse() to obtain the required columns",
      call. = FALSE
    )
  }
  if (any(is.na(raw))) {
    warning("Missing values in '", deparse(substitute(se)), "'")
  }

  if (!is.null(control)) {
    # Show error if control input is not valid
    assertthat::assert_that(
      is.character(control),
      length(control) == 1
    )
    if (!control %in% unique(col_data$condition)) {
      stop("run test_diff() with a valid control.\nValid controls are: '",
        paste0(unique(col_data$condition), collapse = "', '"), "'",
        call. = FALSE
      )
    }
  }

  # variables in formula
  variables <- terms.formula(design_formula) %>%
    attr("variables") %>%
    as.character() %>%
    .[-1]

  # Throw error if variables are not col_data columns
  if (any(!variables %in% colnames(col_data))) {
    stop("run make_diff() with an appropriate 'design_formula'")
  }
  if (variables[1] != "condition") {
    stop("first factor of 'design_formula' should be 'condition'")
  }

  # Obtain variable factors
  for (var in variables) {
    temp <- factor(col_data[[var]])
    assign(var, temp)
  }

  # Make an appropriate design matrix
  design <- model.matrix(design_formula, data = environment())
  colnames(design) <- gsub("condition", "", colnames(design))

  # Generate contrasts to be tested
  # Either make all possible combinations ("all"),
  # only the contrasts versus the control sample ("control") or
  # use manual contrasts
  conditions <- as.character(unique(condition))
  if (type == "all") {
    # All possible combinations
    cntrst <- apply(utils::combn(conditions, 2), 2, paste, collapse = " - ")

    if (!is.null(control)) {
      # Make sure that contrast containing
      # the control sample have the control as denominator
      flip <- grep(paste("^", control, sep = ""), cntrst)
      if (length(flip) >= 1) {
        cntrst[flip] <- cntrst[flip] %>%
          gsub(paste(control, "- ", sep = " "), "", .) %>%
          paste(" - ", control, sep = "")
      }
    }
  }
  if (type == "control") {
    # Throw error if no control argument is present
    if (is.null(control)) {
      stop("run test_diff(type = 'control') with a 'control' argument")
    }

    # Make contrasts
    cntrst <- paste(conditions[!conditions %in% control],
      control,
      sep = " - "
    )
  }
  if (type == "manual") {
    # Throw error if no test argument is present
    if (is.null(test)) {
      stop("run test_diff(type = 'manual') with a 'test' argument")
    }
    assertthat::assert_that(is.character(test))

    if (any(!unlist(strsplit(test, "_vs_")) %in% conditions)) {
      stop("run test_diff() with valid contrasts in 'test'",
        ".\nValid contrasts should contain combinations of: '",
        paste0(conditions, collapse = "', '"),
        "', for example '", paste0(conditions[1], "_vs_", conditions[2]),
        "'.",
        call. = FALSE
      )
    }

    cntrst <- gsub("_vs_", " - ", test)
  }
  # Print tested contrasts
  message(
    "Tested contrasts: ",
    paste(gsub(" - ", "_vs_", cntrst), collapse = ", ")
  )

  # Test for differential expression by empirical Bayes moderation
  # of a linear model on the predefined contrasts
  fit <- limma::lmFit(raw, design = design)
  made_contrasts <- limma::makeContrasts(contrasts = cntrst, levels = design)
  contrast_fit <- limma::contrasts.fit(fit, made_contrasts)

  if (any(is.na(raw))) {
    for (i in cntrst) {
      covariates <- strsplit(i, " - ") %>% unlist()
      single_contrast <- limma::makeContrasts(contrasts = i, levels = design[, covariates])
      single_contrast_fit <- limma::contrasts.fit(fit[, covariates], single_contrast)
      contrast_fit$coefficients[, i] <- single_contrast_fit$coefficients[, 1]
      contrast_fit$stdev.unscaled[, i] <- single_contrast_fit$stdev.unscaled[, 1]
    }
  }

  eB_fit <- limma::eBayes(contrast_fit)

  # function to retrieve the results of
  # the differential expression test using 'fdrtool'
  retrieve_fun <- function(comp, fit = eB_fit) {
    res <- limma::topTable(fit,
      sort.by = "t", adjust.method = "BH", coef = comp,
      number = Inf, confint = TRUE
    )
    # res <- res[!is.na(res$t),]
    # fdr_res <- fdrtool::fdrtool(res$t, plot = FALSE, verbose = FALSE)
    # res$qval <- res$adj.P.Value
    # res$lfdr <- fdr_res$lfdr
    res$comparison <- rep(comp, dim(res)[1])
    res <- tibble::rownames_to_column(res)
    return(res)
  }

  # limma_res<- topTable(eB_fit, sort.by = 'B', adjust.method="BH", coef = cntrst, number = Inf, confint = T )
  # limma_res$comparison <- rep(cntrst, dim(limma_res)[1])
  # limma_res <- rownames_to_column(limma_res)
  # Retrieve the differential expression test results
  limma_res <- purrr::map_df(cntrst, retrieve_fun)

  # Select the logFC, CI and qval variables
  table <- limma_res %>%
    dplyr::select(rowname, logFC, CI.L, CI.R, P.Value, adj.P.Val, comparison) %>%
    dplyr::mutate(comparison = gsub(" - ", "_vs_", comparison)) %>%
    tidyr::gather(variable, value, -c(rowname, comparison)) %>%
    dplyr::mutate(variable = dplyr::recode(variable, logFC = "diff", P.Value = "p.val", adj.P.Val = "p.adj")) %>%
    tidyr::unite(temp, comparison, variable) %>%
    tidyr::spread(temp, value)

  # avoid wrong order of similar comparison names
  comp_list <- sort(gsub(" - ", "_vs_", cntrst))
  ordered_colNames <- c(
    "rowname",
    lapply(comp_list, function(x) {
      colnames(table)[grep(paste0(x, "(_CI.L|_CI.R|_diff|_p.adj|_p.val)"), colnames(table))]
    }) %>% unlist()
  )

  table <- table %>% dplyr::select(dplyr::all_of(ordered_colNames))

  SummarizedExperiment::rowData(se) <- merge(SummarizedExperiment::rowData(se), table,
    by.x = "name", by.y = "rowname", all.x = TRUE
  )
  return(se)
  # return(table)
}

#' @export
get_results_proteins <- function(dep) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"))

  row_data <- SummarizedExperiment::rowData(dep)
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

  # Obtain average protein-centered enrichment values per condition
  row_data$mean <- rowMeans(SummarizedExperiment::assay(dep), na.rm = TRUE)
  centered <- SummarizedExperiment::assay(dep) - row_data$mean
  centered <- data.frame(centered) %>%
    tibble::rownames_to_column() %>%
    tidyr::gather("ID", "val", -rowname) %>%
    dplyr::left_join(data.frame(SummarizedExperiment::colData(dep)), by = "ID")
  centered <- dplyr::group_by(centered, rowname, condition) %>%
    dplyr::summarize(val = mean(val, na.rm = TRUE)) %>%
    dplyr::mutate(val = signif(val, digits = 3)) %>%
    tidyr::spread(condition, val)
  colnames(centered)[2:ncol(centered)] <-
    paste(colnames(centered)[2:ncol(centered)], "_centered", sep = "")

  # Obtain average enrichments of conditions versus the control condition
  ratio <- as.data.frame(row_data) %>%
    # tibble::column_to_rownames("name") %>%
    dplyr::select(dplyr::ends_with("diff")) %>%
    signif(digits = 3) %>%
    tibble::rownames_to_column()
  colnames(ratio)[2:ncol(ratio)] <-
    gsub("_diff", "_log2 fold change", colnames(ratio)[2:ncol(ratio)])
  # df <- left_join(ratio, centered, by = "rowname")

  # Select the adjusted p-values and significance columns
  pval <- as.data.frame(row_data) %>%
    # tibble::column_to_rownames("name") %>%
    dplyr::select(
      dplyr::ends_with("p.val"),
      dplyr::ends_with("p.adj"),
      dplyr::ends_with("significant")
    ) %>%
    tibble::rownames_to_column()
  pval[, grep("p.adj", colnames(pval))] <-
    pval[, grep("p.adj", colnames(pval))] %>%
    signif(digits = 3)
  pval[, grep("p.val", colnames(pval))] <-
    pval[, grep("p.val", colnames(pval))] %>%
    signif(digits = 3)

  # Join into a results table
  ids <- as.data.frame(row_data) %>% dplyr::select(name, ID)
  table <- dplyr::left_join(ids, ratio, by = c("name" = "rowname"))
  table <- dplyr::left_join(table, pval, by = c("name" = "rowname"))
  # table <- dplyr::left_join(table, centered, by = c("name" = "rowname")) %>%
  #   dplyr::arrange(desc(significant))
  table <- as.data.frame(row_data)[, colnames(row_data) %in% c("name", "imputed", "num_NAs", "Protein.names")] %>%
    dplyr::left_join(table, ., by = "name")
  table <- table %>% dplyr::arrange(dplyr::desc(significant))
  colnames(table)[1] <- c("Gene Name")
  colnames(table)[2] <- c("Protein IDs")
  table <- table %>% dplyr::relocate(Protein.names, .after = dplyr::last_col())
  table <- table %>% dplyr::select(grep("[^Protein.names]", colnames(table)), "Protein.names")
  # table$Gene_name<-table$name
  table
}



#######################################################
## Plot Enrichment Results
#######################################################

#' @export
plot_enrichment <- function(gsea_results, number = 10, alpha = 0.05,
                            contrasts = NULL, databases = NULL,
                            nrow = 1, term_size = 8) {
  assertthat::assert_that(
    is.data.frame(gsea_results),
    is.numeric(number),
    length(number) == 1,
    is.numeric(alpha),
    length(alpha) == 1,
    is.numeric(term_size),
    length(term_size) == 1,
    is.numeric(nrow),
    length(nrow) == 1
  )

  # Check gsea_results object
  if (any(!c(
    "Term", "var",
    "contrast", "Adjusted.P.value"
  )
  %in% colnames(gsea_results))) {
    stop("'", deparse(substitute(gsea_results)),
      "' does not contain the required columns",
      "\nMake sure that HGNC gene symbols are present",
      "\n in your 'Gene Names' column of Results table",
      call. = FALSE
    )
  }

  no_enrichment_text <- paste(
    "Enrichment could not be performed.\n",
    "\nDownload enrichment result table for more details. \n"
  )

  if (!is.null(contrasts)) {
    assertthat::assert_that(is.character(contrasts))

    valid_contrasts <- unique(gsea_results$contrast)

    if (!all(contrasts %in% valid_contrasts)) {
      # valid_cntrsts_msg <- paste0("Valid contrasts are: '",
      #                             paste0(valid_contrasts, collapse = "', '"),
      #                             "'")
      # stop("Not a valid contrast, please run `plot_gsea()`",
      #      "with a valid contrast as argument\n",
      #      valid_cntrsts_msg,
      #      call. = FALSE)
      return(
        ggplot2::ggplot() +
          ggplot2::annotate("text", x = 4, y = 25, size = 8, label = no_enrichment_text) +
          ggplot2::theme_void()
      )
    }
    if (!any(contrasts %in% valid_contrasts)) {
      contrasts <- contrasts[contrasts %in% valid_contrasts]
      message(
        "Not all contrasts found",
        "\nPlotting the following contrasts: '",
        paste0(contrasts, collapse = "', '"), "'"
      )
    }

    gsea_results <- dplyr::filter(gsea_results, contrast %in% contrasts)
  }
  if (!is.null(databases)) {
    assertthat::assert_that(is.character(databases))

    valid_databases <- unique(gsea_results$var)

    if (all(!databases %in% valid_databases)) {
      valid_cntrsts_msg <- paste0(
        "Valid databases are: '",
        paste0(valid_databases, collapse = "', '"),
        "'"
      )
      stop("Not a valid database, please run `plot_gsea()`",
        "with valid databases as argument\n",
        valid_cntrsts_msg,
        call. = FALSE
      )
    }
    if (any(!databases %in% valid_databases)) {
      databases <- databases[databases %in% valid_databases]
      message(
        "Not all databases found",
        "\nPlotting the following databases: '",
        paste0(databases, collapse = "', '"), "'"
      )
    }

    gsea_results <- dplyr::filter(gsea_results, var %in% databases)
  }

  # Get top enriched gene sets
  terms <- gsea_results %>%
    dplyr::group_by(contrast, var) %>%
    dplyr::filter(Adjusted.P.value <= alpha) %>%
    dplyr::arrange(Adjusted.P.value) %>%
    dplyr::slice(seq_len(number)) %>%
    .$Term
  subset <- gsea_results %>%
    dplyr::filter(Term %in% terms) %>%
    dplyr::arrange(var, Adjusted.P.value)

  subset$Term <- readr::parse_factor(subset$Term, levels = unique(subset$Term))
  subset$var <- readr::parse_factor(subset$var, levels = unique(subset$var))

  # Plot top enriched gene sets
  if (nrow(subset) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 4, y = 25, size = 8, label = no_enrichment_text) +
      ggplot2::theme_void()
  } else {
    p <- ggplot2::ggplot(subset, ggplot2::aes(Term,
      y = -log10(`Adjusted.P.value`)
    )) +
      ggplot2::geom_col(ggplot2::aes(fill = log_odds)) +
      ggplot2::facet_wrap(~contrast, nrow = nrow) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        y = "-Log10 adjusted p-value",
        fill = "Log2 odds ratio (vs. current background)"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "top",
        legend.text = ggplot2::element_text(size = 9)
      ) +
      ggplot2::scale_fill_distiller(palette = "Spectral")
  }
  p
}

#### ==== get prefix function

get_prefix <- function(words) {
  # Show error if input is not the required class
  assertthat::assert_that(is.character(words))

  # Show error if 'words' contains 1 or less elements
  if (length(words) <= 1) {
    stop("'words' should contain more than one element")
  }
  # Show error if 'words' contains NA
  if (any(is.na(words))) {
    stop("'words' contains NAs")
  }

  # Truncate words to smallest name
  minlen <- min(nchar(words))
  truncated <- substr(words, 1, minlen)

  # Show error if one of the elements is shorter than one character
  if (minlen < 1) {
    stop("At least one of the elements is too short")
  }

  # Get identifical characters
  mat <- data.frame(strsplit(truncated, ""), stringsAsFactors = FALSE)
  identical <- apply(mat, 1, function(x) length(unique(x)) == 1)

  # Obtain the longest common prefix
  prefix <- as.logical(cumprod(identical))
  paste(mat[prefix, 1], collapse = "")
}

#### ===== delete prefix function

delete_prefix <- function(words) {
  # Get prefix
  prefix <- get_prefix(words)
  # Delete prefix from words
  gsub(paste0("^", prefix), "", words)
}

### Filter missing values use different threshold per conditions/groups
threshold_detect <- function(sample_rep) {
  valid_keep <- trunc(sample_rep / 2) + 1
  threshold <- sample_rep - valid_keep
  threshold
}

keep_function <- function(se) {
  # Show error if inputs are not the required classes

  # Show error if inputs do not contain required columns
  if (any(!c("name", "ID") %in% colnames(SummarizedExperiment::rowData(se, use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '",
      deparse(substitute(se)),
      "'\nRun make_unique() and make_se() to obtain the required columns",
      call. = FALSE
    )
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(SummarizedExperiment::colData(se)))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
      deparse(substitute(se)),
      "'\nRun make_se() or make_se_parse() to obtain the required columns",
      call. = FALSE
    )
  }

  # Make assay values binary (1 = valid value)
  bin_data <- SummarizedExperiment::assay(se)

  # new, removed ref columns from dataset
  bin_data <- bin_data %>% data.frame()
  idx <- is.na(bin_data) # idx <- is.na(SummarizedExperiment::assay(se))
  bin_data[!idx] <- 1
  bin_data[idx] <- 0

  # Filter se on the maximum allowed number of
  # missing values per condition (defined by thr)
  keep <- bin_data %>%
    data.frame() %>%
    tibble::rownames_to_column() %>%
    tidyr::gather("ID", "value", -rowname) %>%
    dplyr::left_join(data.frame(SummarizedExperiment::colData(se)), by = "ID") %>%
    dplyr::group_by(rowname, condition) %>%
    dplyr::summarize(miss_val = dplyr::n() - sum(value))
  keep
}

#' @export
filter_missval_new <- function(se, one_condition, exp_df) {
  threshold <- exp_df$thr[exp_df$condition == one_condition] %>% unlist()
  keep <- keep_function(se)
  keep1 <- keep %>%
    dplyr::filter(condition == one_condition) %>%
    dplyr::filter(miss_val <= threshold)
  keep1 <- keep1 %>%
    tidyr::spread(condition, miss_val)
  se_fltrd <- se[keep1$rowname, ]
  se_fltrd
}
