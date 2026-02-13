#' Generate multi-omic differential networks
#'
#' Generate a multi-layer differential network with interaction p values
#'
#' @param assayData a matrix (or list of matrices for multi-omic analysis) 
#' containing normalised assay data. Sample IDs must be in columns and probe 
#' IDs (genes, proteins...) in rows. The matrix must be numeric.
#' For multi-omic analysis, it is highly recommended to use a named list of
#' matrices. 
#' If unnamed, sequential names (assayData1, assayData2, etc.) will 
#' be assigned to identify each matrix.
#' @param metadata a named vector, matrix, or data.frame containing sample 
#' annotations or categories. If matrix or data.frame, each row should
#' correspond to a sample, with columns representing different sample 
#' characteristics (e.g., treatment group, condition, time point). The colname 
#' of the sample characteristic to be used for differential analysis must be 
#' specified in `category_variable`. Rownames must match the sample IDs used in
#' assayData. 
#' If named vector, each element must correspond to a sample characteristic to 
#' be used for differential analysis, and names must match sample IDs used in
#' the colnames of `assayData`.
#' Continuous variables are not allowed.  
#' @param category_variable when metadata is a matrix or data.frame this is the 
#' column name of `metadata` that contains the sample annotations to be used for
#' differential analysis
#' @param regression_method whether to use robust linear modelling to calculate
#' link p values. Options are 'lm' (default) or 'rlm'. The lm implementation 
#' is faster and lighter.  
#' @param mixedModel logical. If TRUE, use mixed models with random intercept
#' for id variable. Default FALSE.
#' @param id_variable character. Name of the column in metadata containing the
#' grouping variable for random effects (e.g., "slide_id", "batch"). Only used
#' when mixedModel = TRUE. Default NULL.
#' @param lmer_ctrl list of control parameters for lmer mixed models.
#' Only used when mixedModel = TRUE. See ?lme4::lmerControl for details.
#' @param category_subset optional character vector indicating a subset of 
#' categories from the category variable. If not specified, all categories in
#' `category_variable` will be used.
#' @param network network of biological interactions provided by the user. The
#' network must be provided in the form of a table of class data.frame with only 
#' two columns named "from" and "to".
#' If NULL (default) a network of 10,537 molecular interactions obtained from
#' KEGG, mirTARbase, miRecords and transmiR will be used.
#' This has been obtained via the `exportgraph` function of the MITHrIL tool
#' (Alaimo et al., 2016).
#' @param  percentile_vector a numeric vector specifying the percentiles to be
#' used in the percolation analysis. By default, it is defined as 
#' `seq(0.35, 0.98, by = 0.05)`, which generates a sequence of percentiles
#' starting at 0.35, meaning that targets (genes/proteins...) whose expression 
#' value is under the 35th percentile of the whole matrix will be excluded.
#' This threshold can be modified by specifying a different starting point for 
#' `seq`. For a more granular percolation analysis an higher optimisation of 
#' the algorithm, `by = 0.05` can be modified in favour of lower values, but 
#' this will increase the computational time.  
#' @param padj_method a character string indicating the p values correction 
#' method for multiple test adjustment. It can be either one of the methods 
#' provided by the `p.adjust` function from `stats` (bonferroni, BH, hochberg, 
#' etc.) or "q.value" for Storey's q values, or "none" for unadjusted p values.
#' When using "q.value" the `qvalue` package must be installed first.
#' @param show_progressBar logical. Whether to display a progress bar during 
#' execution. Default is TRUE.
#' @param verbose logical. Whether to print detailed output messages during 
#' processing. Default is TRUE
#' @param cores number of cores to use for parallelisation.
#' @return a `deggs` object containing differential networks incorporating
#' p values or adjusted p values for each link.
#' @examples
#' # Single omic analysis:
#' data("synthetic_metadata")
#' data("synthetic_rnaseqData")
#' deggs_object <- get_diffNetworks(assayData = synthetic_rnaseqData,
#'                                  metadata = synthetic_metadata,
#'                                  category_variable = "response",
#'                                  regression_method = "lm",
#'                                  padj_method = "bonferroni",
#'                                  verbose = FALSE,
#'                                  show_progressBar = FALSE,
#'                                  cores = 1)
#'  
#' # multi-omic analysis: 
#' data("synthetic_metadata")
#' data("synthetic_rnaseqData")
#' data("synthetic_proteomicData")
#' data("synthetic_OlinkData")
#' assayData_list <- list("RNAseq" = synthetic_rnaseqData,
#'                        "Proteomics" = synthetic_proteomicData,
#'                        "Olink" = synthetic_OlinkData)
#' deggs_object <- get_diffNetworks(assayData = assayData_list,
#'                                  metadata = synthetic_metadata,
#'                                  category_variable = "response",
#'                                  regression_method = "lm",
#'                                  padj_method = "bonferroni",
#'                                  verbose = FALSE,
#'                                  show_progressBar = FALSE,
#'                                  cores = 1)
#'                                  
#' # to use only certain categories for comparison: 
#' # let's randomly add another level of response to the example metadata
#' synthetic_metadata$response <- as.character(synthetic_metadata$response)
#' indices <- sample(1:nrow(synthetic_metadata), 20, replace = FALSE) 
#' synthetic_metadata$response[indices] <- "Moderate response"
#' deggs_object <- get_diffNetworks(assayData = assayData_list,
#'                                  metadata = synthetic_metadata,
#'                                  category_variable = "response",
#'                                  category_subset = c("Responder", 
#'                                                      "Non_responder"),
#'                                  regression_method = "lm",
#'                                  verbose = FALSE,
#'                                  show_progressBar = FALSE,
#'                                  cores = 1)
#'                                  
#' # to be more generous on the targets to be excluded, and lower the expression 
#' # level threshold to the 25th percentile (or lower): 
#' deggs_object <- get_diffNetworks(assayData = assayData_list,
#'                                  metadata = synthetic_metadata,
#'                                  category_variable = "response",
#'                                  category_subset = c("Responder", 
#'                                                      "Non_responder"),
#'                                  regression_method = "lm",
#'                                  percentile_vector = seq(0.25, 0.98, by = 0.05),
#'                                  verbose = FALSE,
#'                                  show_progressBar = FALSE,
#'                                  cores = 1)                                  
#' @export
get_diffNetworks <- function(assayData,
                             metadata,
                             category_variable = NULL,
                             regression_method = 'lm',
                             mixedModel = FALSE,
                             id_variable = NULL,
                             lmer_ctrl = NULL,
                             category_subset = NULL,
                             network = NULL,
                             percentile_vector = seq(0.35, 0.98, by = 0.05),
                             padj_method = "bonferroni",
                             show_progressBar = TRUE,
                             verbose = TRUE, 
                             cores = parallel::detectCores() / 2) {
  
  regression_method <- match.arg(regression_method, choices = c('lm', 'rlm'))
  
  # Check assayData validity (matrix or list)
  if (!(is.list(assayData) || is.matrix(assayData))) {
    stop("assayData must be a matrix or a list of matrices")
  }
  
  # Convert single matrix to list
  if (is.matrix(assayData)) {
    assayDataList <- list(assayData)
  } else {
    # If it's a list, check that all elements are matrices
    if (!all(vapply(assayData, is.matrix, logical(1)))) {
      non_matrix <- which(!vapply(assayData, is.matrix, logical(1)))
      stop("All elements of assayData list must be matrices. ",
           "Elements ", paste(non_matrix, collapse = ", "), " are not matrices.")
    }
    assayDataList <- assayData
  }
  
  # check that all matrices are numeric
  check_numeric <- vapply(assayDataList, function(x) is.numeric(x), logical(1))
  if (!all(check_numeric)) {
    non_numeric <- which(!check_numeric)
    stop("All matrices in assayData must be numeric. ",
         "Elements ", paste(non_numeric, collapse = ", "), " are not numeric.")
  }
  
  # Check category_variable validity
  if (!is.null(category_variable)) {
    # If metadata is a matrix/data.frame, category_variable must be a column name
    if ((is.matrix(metadata) || is.data.frame(metadata)) && 
        !category_variable %in% colnames(metadata)) {
      stop("category_variable must be %in% colnames(metadata). If metadata is a 
       named vector category_variable must be NULL.")
    }
  } else if (is.matrix(metadata) || is.data.frame(metadata)) {
    # If metadata is a matrix/data.frame, category_variable cannot be NULL
    stop("metadata is a matrix or dataframe but category_variable is NULL.
     Please specify category_variable to select the correct column 
     from metadata.")
  }
  
  # Check id_variable for mixed models 
  if (mixedModel && is.null(id_variable)) {
    stop("id_variable must be provided when mixedModel = TRUE")
  }
  
  # Check required packages for mixed models
  if (mixedModel) {
    if (regression_method == 'lm') {
      if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("Package 'lme4' is required for mixed models. Please install it.")
      }
    } else {
      if (!requireNamespace("robustlmm", quietly = TRUE)) {
        stop("Package 'robustlmm' is required for robust mixed models. Please install it.")
      }
    }
  }
  
  # set controls for mixed models 
  if (mixedModel && regression_method == "lm") {
    # Control settings for lmer
    lmer_ctrl <- lme4::lmerControl(
      optimizer = "bobyqa",  # Fast and generally reliable
      check.conv.singular = "ignore",  # I handle this separately
    )
  }
  
  if (is.null(names(assayDataList))) (
    names(assayDataList) <- paste0("assayData", seq_along(assayDataList))
  )
  
  tidy_data <- tidy_metadata(
    category_subset = category_subset,
    metadata = metadata,
    category_variable = category_variable,
    id_variable = id_variable,
    verbose = verbose
  )  
  
  if (is.list(tidy_data)) {
    metadata <- tidy_data[["metadata"]] 
    id_variable <- tidy_data[["id_vector"]] 
  } else {
    metadata <- tidy_data  
  }
  # metadata (and, when present id_variable) will be a named vector from now on
  
  diffNetworks_list <- lapply(seq_along(assayDataList), function(i){
    assayDataName <- names(assayDataList)[i]
    assayData <- assayDataList[[assayDataName]]
    
    diffNetworks <- 
      tryCatch({
        get_diffNetworks_singleOmic(assayData = assayData,
                                    assayDataName = assayDataName, 
                                    metadata = metadata,
                                    regression_method = regression_method,
                                    mixedModel = mixedModel,
                                    id_variable = id_variable,
                                    lmer_ctrl = lmer_ctrl, 
                                    network = network,
                                    percentile_vector = percentile_vector,
                                    padj_method = padj_method,
                                    show_progressBar = show_progressBar,
                                    verbose = verbose, 
                                    cores = cores)
      }, error = function(e) { message(conditionMessage(e)) })
  })
  names(diffNetworks_list) <- names(assayDataList)
  
  degg <- list(
    diffNetworks = diffNetworks_list,
    assayData = assayDataList,
    metadata = metadata,
    regression_method = regression_method,
    mixedModel = mixedModel,
    id_variable = id_variable,
    lmer_ctrl = lmer_ctrl,
    category_subset = category_subset,
    padj_method = padj_method
  )
  class(degg) <- "deggs"
  return(degg)
}

#' Generate differential networks for single omic analysis 
#' 
#' @inheritParams get_diffNetworks
#' @param assayDataName name of the assayData, to identify which omic is. 
#' @return a list of differential networks, one per category
get_diffNetworks_singleOmic <- function(assayData,
                                        assayDataName,
                                        metadata,
                                        regression_method,
                                        mixedModel,
                                        id_variable,
                                        lmer_ctrl,
                                        network,
                                        percentile_vector, 
                                        padj_method,
                                        show_progressBar,
                                        verbose, 
                                        cores) {
  
  if (verbose) (
    message("Generating ", assayDataName, " network layer...")
  )
  
  if (!(is.matrix(assayData) || is.data.frame(assayData))) (
    stop(assayDataName, " is not a matrix or data.frame")
  )
  
  if(!any(colnames(assayData) %in% names(metadata))) (
    stop("Sample IDs in ", assayDataName, " don't match the IDs in metadata 
          (the metadata names/rownames must match the sample IDs used as column 
          names of ", assayDataName, ".)")
  )
  
  # align metadata and assayData
  if (length(intersect(colnames(assayData), names(metadata))) == 0) {
    stop("Sample IDs in metadata and ", assayDataName, " don't match.")
  }
  
  # remove assayData samples that don't exist in metadata (because they
  # were missing or because they were filtered out due to the category_subset)
  assayData <- assayData[, which(colnames(assayData) %in% names(metadata)),
                         drop = FALSE]
  # remove cols with 0 or 1 non-NA values (we need at least two points)
  assayData <- assayData[, colSums(is.na(assayData)) < (nrow(assayData) - 1)]
  
  # message on metadata sample IDs that don't exist in assayData
  missing_metadataSamples <- names(metadata)[which(!(names(metadata) %in% 
                                                       colnames(assayData)))]
  
  if (length(missing_metadataSamples) != 0) {
    missing_samples_formatted <- paste(missing_metadataSamples, collapse = ", ")
    message("The following samples IDs are missing in ", assayDataName, 
            ":\n", missing_samples_formatted)
  }
  
  # align 
  metadata <- metadata[colnames(assayData)]
  
  if(length(unique(metadata)) == 1) (
    stop("All sample IDs in ", assayDataName, " belong to one category. 
           No differential analysis is possible.")
  )
  
  if (!is.null(network)) (
    # user provided network
    network_to_use <- network
  ) else (
    network_to_use <- metapathway_gene_symbols
  )
  
  edges <- network_to_use[network_to_use$from %in% rownames(assayData) & 
                            network_to_use$to %in% rownames(assayData), ] 
  
  if (nrow(edges) == 0) (
    stop("Rownames of ", assayDataName, " don't match any link of the biological
    network. Check if names need conversion to gene symbols or entrez IDs. 
    Alternatively, provide more data.")
  )
  
  # Remove duplicate edges
  edges$edge_ID <- paste(edges$from, edges$to)
  edges <- edges[!duplicated(edges$edge_ID), ]
  edges$edge_ID <- NULL

  nodes <- c(edges[, 1], edges[, 2]) %>%
    unique()
  
  assayData <- assayData[rownames(assayData) %in% nodes, ]
  
  categories <- levels(metadata)
  
  # pairwise contrasts are necessary for more than 2 categories  
  contrasts <- utils::combn(categories, m = 2) %>%
    as.data.frame()
  
  # create categories (duplicating assay data for each category)
  category_median_list <- lapply(categories, function(one_category) {
    # Get sample IDs belonging to the current category
    sample_ids_in_category <- names(metadata)[metadata == one_category]
    
    # Select samples from assayData using these IDs
    category_vec <- assayData[, sample_ids_in_category, drop = FALSE]
    
    # Calculate median expression across samples
    category_vec <- apply(category_vec, 1, stats::median, na.rm = TRUE)
    rownames(category_vec) <- NULL
    return(category_vec)
  })
  names(category_median_list) <- categories
  
  # Percolation analysis: 
  # calculate interaction p values per percentile and select the percentile 
  # network with the highest number of significant interactions 
  sig_edges_count <- 0

  if (Sys.info()["sysname"] == "Windows") {
    cl <- parallel::makeCluster(cores)
    
    parallel::clusterExport(cl, c(
      "percentile_vector", "category_median_list", "contrasts",
      "regression_method", "edges", "categories", "calc_pvalues_percentile",
      "assayData", "metadata",
      "sig_edges_count", "padj_method"
    ), envir = environment())
    
    if (show_progressBar) (
      pvalues_list <- pbapply::pblapply(
        cl = cl, percentile_vector,
        function(percentile) {
          calc_pvalues_percentile(
            assayData = assayData,
            padj_method = padj_method,
            metadata = metadata,
            percentile = percentile,
            category_median_list = category_median_list,
            contrasts = contrasts,
            regression_method = regression_method,
            mixedModel = mixedModel,
            id_variable = id_variable,
            lmer_ctrl = lmer_ctrl,
            edges = edges,
            categories_length = length(categories),
            sig_edges_count = sig_edges_count
          )
        })
    ) else (
      pvalues_list <- parallel::parLapply(
        cl = cl, percentile_vector,
        function(percentile) {
          calc_pvalues_percentile(
            assayData = assayData,
            padj_method = padj_method,
            metadata = metadata,
            percentile = percentile,
            category_median_list = category_median_list,
            contrasts = contrasts,
            regression_method = regression_method,
            mixedModel = mixedModel,
            id_variable = id_variable,
            lmer_ctrl = lmer_ctrl,
            edges = edges,
            categories_length = length(categories),
            sig_edges_count = sig_edges_count
          )
        })
    )
    parallel::stopCluster(cl)
  } else {
    # Parallelisation for mac/linux
    if (show_progressBar) (
      pvalues_list <- pbmcapply::pbmclapply(
        mc.cores = cores, percentile_vector,
        function(percentile) (
          calc_pvalues_percentile(
            assayData = assayData,
            padj_method = padj_method,
            metadata = metadata,
            percentile = percentile,
            category_median_list = category_median_list,
            contrasts = contrasts,
            regression_method = regression_method,
            mixedModel = mixedModel,
            id_variable = id_variable,
            lmer_ctrl = lmer_ctrl,
            edges = edges,
            categories_length = length(categories),
            sig_edges_count = sig_edges_count
          )
        )
      )
    ) else (
      pvalues_list <- parallel::mclapply(
        mc.cores = cores, percentile_vector,
        function(percentile) (
          calc_pvalues_percentile(
            assayData = assayData,
            padj_method = padj_method,
            metadata = metadata,
            percentile = percentile,
            category_median_list = category_median_list,
            contrasts = contrasts,
            regression_method = regression_method,
            mixedModel = mixedModel,
            id_variable = id_variable,
            lmer_ctrl = lmer_ctrl,
            edges = edges,
            categories_length = length(categories),
            sig_edges_count = sig_edges_count
          )
        )
      )
    )
    
  }
  
  # Extract values if warnings are present (for mc.cores = 1)
  if (cores == 1 && !is.null(pvalues_list$value)) {
    if (!is.null(pvalues_list$warning)) {
      warning(pvalues_list$warning$message, call. = FALSE)
    }
    pvalues_list <- pvalues_list$value
  }
  
  names(pvalues_list) <- percentile_vector
  
  # Extract the network filtered at the optimal threshold percentile
  # (the one with the highest number of statistically significant interactions)
  sig_pvalues <- vapply(pvalues_list, function(networks) {
    if (!inherits(networks, "try-error") && !is.null(networks)) {
      return(networks$sig_pvalues_count)
    } else {
      return(NA_real_)
    }
  }, FUN.VALUE = numeric(1))
  
  sig_pvalues <- sig_pvalues[!is.na(sig_pvalues)]
  
  if (length(sig_pvalues) == 0) {
    stop("No significant differences detected across categories.")
  }
  
  best_percentile <- sig_pvalues[which.max(sig_pvalues)]
  if (verbose) (
    message("Percolation analysis: genes whose expression is below the ",
            as.numeric(names(best_percentile)) * 100,
            "th percentile are removed from networks.")
  )
  final_networks <- pvalues_list[[as.character(names(best_percentile))]]

  # Report singularity issues if mixed model was used
  if (mixedModel) {
    if (verbose) {
      n_singular <- sum(sapply(final_networks, function(network) {
        if (is.data.frame(network)) sum(network$singular) else 0
      }))
      
      n_total <- sum(sapply(final_networks, function(network) {
        if (is.data.frame(network)) nrow(network) else 0
      }))
      
      perc_singular <- round(100 * n_singular / n_total, 1)
      
      if (n_singular > 0) {
        message(sprintf(
          "Singularity detected in %d out of %d model fits (%.1f%%).",
          n_singular, n_total, perc_singular
        ))
        message("Non-mixed models were used as fallback for singular fits.")
        
        if (perc_singular > 50) {
          warning(
            "More than 50% of models were singular. ",
            "Consider using non-mixed models (mixedModel = FALSE) for this analysis."
          )
        }
      }
    }
  }
  return(final_networks)
}


#' Tidying up of metadata and id_variable (when provided). 
#' 
#' Samples belonging to undesidered categories
#' (if specified) will be removed as well as categories with less than five 
#' samples, and NAs.
#' Both metadata and id_variable (when provided) will be returned as named factors. 
#' Used internally by get_diffNetworks_singleOmic.
#'
#' @inheritParams get_diffNetworks
#' @return a tidy named factor vector of sample annotations, 
#' or a list with metadata and id_vector when mixedModel = TRUE
#' @keywords internal
tidy_metadata <- function(category_subset,
                          metadata,
                          category_variable,
                          id_variable,
                          verbose) {
  
  # Convert metadata to a named vector based on input type
  if (is.atomic(metadata) || is.factor(metadata)) {
    if (is.null(names(metadata))) {
      stop("Named vector required: metadata vector must have names.")
    }
    # Already a named vector
    category_vector <- metadata
  } else if (is.matrix(metadata) || is.data.frame(metadata)) {
    if (is.null(category_variable)) {
      stop("category_variable must be specified when metadata is
           a matrix or data.frame.")
    }
    if (!category_variable %in% colnames(metadata)) {
      stop("category_variable '", category_variable, "' not found in metadata.")
    }
    category_vector <- metadata[, category_variable]
    names(category_vector) <- rownames(metadata)
  } else {
    stop("metadata must be a named vector/factor, matrix, or data.frame.")
  }
  
  # Extract id_vector if id_variable is provided
  if (!is.null(id_variable)) {
    if (is.atomic(metadata) || is.factor(metadata)) {
      stop("When id_variable is provided, metadata must be a matrix or data.frame.")
    }
    if (!id_variable %in% colnames(metadata)) {
      stop("id_variable '", id_variable, "' not found in metadata.")
    }
    id_vector <- metadata[, id_variable]
    names(id_vector) <- rownames(metadata)
    
    # Check for NAs in id_vector
    if (any(is.na(id_vector))) {
      stop("id_variable contains NA values. Please remove missing values.")
    }
  }
  
  # Subset selected categories (optional)
  if (!is.null(category_subset)) {
    if (!all(category_subset %in% category_vector)) (
      stop("category_subset values must be contained in metadata")
    )
    keep_samples <- names(category_vector)[category_vector %in% category_subset]
    if (length(keep_samples) == 0) {
      stop("No samples remain after filtering for specified category_subset.")
    }
    category_vector <- category_vector[keep_samples]
    if (!is.null(id_variable)) {
      id_vector <- id_vector[keep_samples]
    }
  }
  
  # Convert to factor
  if (!is.factor(category_vector)) {
    category_vector <- as.factor(category_vector)
    if (verbose) {
      message("category values were converted to factor")
    }
  }
  
  if (!is.null(id_variable)) {
    if (!is.factor(id_vector)) {
      id_vector <- as.factor(id_vector)
      if (verbose) {
        message("id values were converted to factor")
      }
    }
  }
  
  # Remove categories with less than five observations
  tbl <- table(category_vector)
  small_groups <- names(tbl)[tbl < 5]
  
  if (length(small_groups) > 0) {
    for (group in small_groups) {
      message("The ", group, 
              " category did not contain enough samples (less than five 
              observations). ",
              "This category will be removed.")
    }
    keep_samples <- names(category_vector)[category_vector %in%
                                             names(tbl)[tbl >= 5]]
    category_vector <- category_vector[keep_samples]
    category_vector <- droplevels(category_vector)
    
    if (!is.null(id_variable)) {
      id_vector <- id_vector[keep_samples]
      id_vector <- droplevels(id_vector)
    }
    
    if (nlevels(category_vector) < 2) {
      stop("At least two categories with at least 5 observations are required.")
    }
  }
  
  # Remove NAs
  NAs_number <- sum(is.na(category_vector))
  if (NAs_number > 0) {
    keep_samples <- names(category_vector)[!is.na(category_vector)]
    category_vector <- category_vector[keep_samples]
    if (!is.null(id_variable)) {
      id_vector <- id_vector[keep_samples]
    }
    if (verbose) {
      message("category values contained NAs. ",
              NAs_number, " samples were removed.")
    }
  }
  
  # Return
  if (is.null(id_variable)) {
    return(category_vector)
  } else {
    return(list("metadata" = category_vector, "id_vector" = id_vector))
  }
}


#' Compute interaction p values for a single percentile value
#'
#' @inheritParams get_diffNetworks_singleOmic
#' @param categories_length integer number indicating the number of categories
#' @param category_median_list list of category data.frames
#' @param percentile a float number indicating the percentile to use.
#' @param contrasts data.frame containing the categories contrasts in rows
#' @param edges network of biological interactions in the form of a table of
#' class data.frame with two columns: "from" and "to".
#' @param sig_edges_count number of significant edges (p < 0.05)
#' @importFrom magrittr %>%
#' @return The list of float numbers of the significant pvalues for a single
#' percentile
calc_pvalues_percentile <- function(assayData,
                                    metadata,
                                    categories_length,
                                    category_median_list,
                                    padj_method,
                                    percentile,
                                    contrasts,
                                    regression_method,
                                    mixedModel,
                                    id_variable,
                                    lmer_ctrl,
                                    edges,
                                    sig_edges_count) {
  
  # 1st filtering step: remove low expressed genes in each category
  # (genes under the percentile threshold)
  cut_off <- stats::quantile(assayData, prob = percentile,
                             na.rm = TRUE)
  user_message <- paste0("No gene above the threshold (", percentile * 100,
                         "th percentile).")
  
  category_median_list <- lapply(category_median_list, function(category_vec) {
    category_vec <- category_vec[category_vec > cut_off]
    if (length(category_vec) == 0) (
      category_vec <- user_message
    )
    return(category_vec)
  })
  
  # format to string vectors to allow easy detection of overlapping edges
  networks_to_string <- lapply(category_median_list, function(category_vec) {
    if (!is.character(category_vec)) {
      categoryEdges <- edges[edges$from %in% names(category_vec) & 
                               edges$to %in% names(category_vec), ]
      network_to_string <- do.call(paste, categoryEdges)
    } else {
      # assign error message for empty vectors
      network_to_string <- user_message
    }
    return(network_to_string)
  })
  
  # find overlapping edges across categories
  common_links <- lapply(contrasts, function(contrast) {
    return(intersect(
      networks_to_string[[contrast[1]]],
      networks_to_string[[contrast[2]]]
    ))
  }) %>%
    unlist() %>%
    unique()
  
  # remove error message from final list links 
  if (user_message %in% common_links) (
    common_links <- common_links[!(common_links %in% user_message)]
  )

  # 2nd filtering step
  # remove overlapping edges and format back to data.frame
  categories_network_list <- lapply(networks_to_string, function(subnetwork) {
    subnetwork <- subnetwork[!subnetwork %in% common_links]
    if (!(user_message %in% subnetwork)) (
      subnetwork <- data.frame(
        'from' = unlist(lapply(strsplit(subnetwork, " "), `[[`, 1)),
        'to' = unlist(lapply(strsplit(subnetwork, " "), `[[`, 2))
      )
    )
    return(subnetwork)
  })
  
  # count tot edges left per category
  tot_edges <- lapply(categories_network_list, function(category_network) {
    if (!is.character(category_network)) (
      nrow(category_network)
    )
  }) %>% 
    unlist() %>% 
    sum()

  # calculate interaction p values
  pvalues_list <- lapply(categories_network_list, function(category_network) {
    return(calc_pvalues_network(
      category_network = category_network,
      assayData = assayData,
      metadata = metadata,
      regression_method = regression_method,
      categories_length = categories_length,
      padj_method = padj_method,
      mixedModel = mixedModel,
      id_variable = id_variable,
      lmer_ctrl = lmer_ctrl
    ))
  })
  # count significant p values
  sig_var <- ifelse(padj_method == "none", "p.value", "p.adj")
  if (tot_edges > sig_edges_count) {
    # num tot edges greater than previous sig edges count
    p_values_sig_count <- lapply(pvalues_list, function(category_network) {
      if (!is.null(dim(category_network))) {
        num_sig_links <- sum(category_network[, sig_var] < 0.05, na.rm = TRUE)
        # return the number of links with significant p values or adj p
        if (!is.na(num_sig_links)) {
          return(num_sig_links)
        } else {return(0)}
      } else {return(0)}
    }) %>% unlist()
    num_sig_pvalues <- sum(p_values_sig_count) # these are p values or adj p
    
        if (num_sig_pvalues > sig_edges_count) {
      sig_edges_count <<- num_sig_pvalues
        }
    # The <<- operator here is used only for efficiency reasons and does not
    # alter the .GlobalEnv as the variable sig_edges_count is defined in the 
    # parent environment of the get_diffNetworks_singleOmic() parent function 
    # (see
    # https://contributor.r-project.org/cran-cookbook/code_issues.html#writing-to-the-.globalenv).  
    # <<- represents the only way to implement early stopping of the loop 
    # when the total number of edges left in the network for a given 
    # percentile is smaller that the maximum number of significant edges found 
    # in previous percentile's networks. It wouldn't make sense to keep looking 
    # for a bigger number of significant edges in subsequent percentiles. 
    # The update of sig_edges_count will allow skipping unnecessary lm or rlm
    # iterations. 

    pvalues_list <- append(pvalues_list, num_sig_pvalues)
    names(pvalues_list)[length(pvalues_list)] <- "sig_pvalues_count"
    return(pvalues_list)
  } else {
    # previous sig edges count greater than num tot edges, skip all calculations
    return(NULL)
  }
}


#' Calculate the p values for specific category network samples
#'
#' @inheritParams calc_pvalues_percentile
#' @param category_network network table for a specific category
#' @importFrom methods is
#' @importFrom MASS rlm
#' @importFrom sfsmisc f.robftest  
#' @return a list of p values
calc_pvalues_network <- function(assayData,
                                 metadata,
                                 padj_method,
                                 categories_length,
                                 regression_method,
                                 category_network,
                                 mixedModel,
                                 id_variable,
                                 lmer_ctrl) {
  
  if (!is.character(category_network)) {
    if (nrow(category_network) > 0) {
      
      # Track singularity and fallback to non-mixed models
      singular <- rep(FALSE, nrow(category_network))
      # For mixed-effects models, 
      # the singular vector will be updated using the <<- operator in the code
      # below. This is done for efficiency reasons and does not
      # alter the .GlobalEnv as the variable singular is defined in the 
      # parent environment of the calc_pvalues_percentile parent function 
      # (see https://contributor.r-project.org/cran-cookbook/code_issues.html#writing-to-the-.globalenv).  

      p.value <- vapply(seq_len(nrow(category_network)), function(i){
        gene_A <- category_network[i, 1]
        gene_B <- category_network[i, 2]
        
        binary.metadata <- as.numeric(metadata) - 1
        
        # Two categories case 
        if (categories_length == 2) {
          
          # Standard linear model
          if (regression_method == "lm") {
            if (!mixedModel) {
              p_interaction <- fit_lm_interaction(
                assayData[gene_A, ], 
                assayData[gene_B, ], 
                binary.metadata
              )
            } else {
              # Mixed linear model
              lmm_fit <- lme4::lmer(assayData[gene_B, ] ~
                                      assayData[gene_A, ] * 
                                      binary.metadata + (1 | id_variable),
                                    control = lmer_ctrl)
              
              # if singular, fallback to standard linear model (and extract p value)
              if (lme4::isSingular(lmm_fit)) {
                singular[i] <<- TRUE
                p_interaction <- fit_lm_interaction(
                  assayData[gene_A, ], 
                  assayData[gene_B, ], 
                  binary.metadata
                )
              # if not singular, extract p value from the fitted mixed model
              } else {
                coef_summary <- summary(lmm_fit)$coefficients
                if (nrow(coef_summary) >= 4) {
                  t_value <- coef_summary[4, 3]
                  dof <- length(binary.metadata) - length(lme4::fixef(lmm_fit))
                  p_interaction <- 2 * (1 - stats::pt(abs(t_value), df = dof))
                } else {
                  p_interaction <- NA_real_
                }
              }
            }
          }
          
          # Robust linear model
          if (regression_method == "rlm") {
            if (!mixedModel) {
              robustfit <- rlm(assayData[gene_B, ] ~
                                 assayData[gene_A, ] * binary.metadata)
              
              p_interaction <- try(
                f.robftest(robustfit, var = 3)$p.value, silent = TRUE
              )
              if (inherits(p_interaction, "try-error")) {
                p_interaction <- NA_real_
              }
            } else {
              # Robust mixed linear model
              rlmm_fit <- suppressMessages(
                robustlmm::rlmer(assayData[gene_B, ] ~
                                   assayData[gene_A, ] * 
                                   binary.metadata + (1 | id_variable))
              )
              
              # if singular, fallback to standard linear model (and extract p value)
              if (lme4::isSingular(rlmm_fit)) {
                singular[i] <<- TRUE
                robustfit <- rlm(assayData[gene_B, ] ~
                                   assayData[gene_A, ] * binary.metadata)
                
                p_interaction <- try(
                  f.robftest(robustfit, var = 3)$p.value, silent = TRUE
                )
                if (inherits(p_interaction, "try-error")) {
                  p_interaction <- NA_real_
                }
              # if not singular, extract p value from the fitted mixed model
              } else {
                coef_summary <- summary(rlmm_fit)$coefficients
                if (nrow(coef_summary) >= 4) {
                  t_value <- coef_summary[4, 3]
                  dof <- length(binary.metadata) - length(lme4::fixef(rlmm_fit))
                  p_interaction <- 2 * (1 - stats::pt(abs(t_value), df = dof))
                } else {
                  p_interaction <- NA_real_
                }
              }
            }
          }
        }
        
        # Three or more categories case 
        if (categories_length >= 3) {
          if (!mixedModel) {
            # Standard ANOVA
            res_aov <- stats::aov(assayData[gene_B, ] ~
                                    assayData[gene_A, ] * metadata)
            p_interaction <- summary(res_aov)[[1]][["Pr(>F)"]][3]
          } else {
            # Mixed model ANOVA
            if (!requireNamespace("car", quietly = TRUE)) {
              stop("Package 'car' is required for mixed models and >2 categories. Please install it.")
            }
            lmm_fit <- suppressMessages(
              lme4::lmer(assayData[gene_B, ] ~
                           assayData[gene_A, ] *
                           metadata + (1 | id_variable),
                         control = lmer_ctrl)
            )
            
            # if singular, fallback to standard linear model (and extract p value)
            if (lme4::isSingular(lmm_fit)) {
              singular[i] <<- TRUE
              res_aov <- stats::aov(assayData[gene_B, ] ~
                                      assayData[gene_A, ] * metadata)
              p_interaction <- summary(res_aov)[[1]][["Pr(>F)"]][3]
              # if not singular, extract p value from the fitted mixed model
            } else {
              anova_table <- car::Anova(lmm_fit, type = "III")
              if (nrow(anova_table) >= 4) {
                p_interaction <- anova_table[4, 3]
              } else {
                p_interaction <- NA_real_
              }
            }
          }
        }
        
        return(p_interaction)
      }, FUN.VALUE = numeric(1))
    } else {
      p.value <- "No specific links for this category."
    }
  } else {
    p.value <- "No specific links for this category."
  }
  
  if (is.numeric(p.value)) {
    category_network <- cbind(category_network, p.value, singular)
    
    if (padj_method == "q.value") {
      p.adj <- try(qvalue::qvalue(category_network[, "p.value"])$qvalues,
                   silent = TRUE)
      if (is(p.adj, "try-error")) {
        if (nrow(category_network) > 1) {
          p.adj <- qvalue::qvalue(p = category_network[, "p.value"],
                                  pi0 = 1)$qvalues
        } else {
          p.adj <- NA
        }
      }
      category_network$p.adj <- p.adj
    }
    
    if (!(padj_method %in% c("q.value", "none"))) {
      p.adj <- try(stats::p.adjust(category_network[, "p.value"],
                                   method = padj_method))
      if (is(p.adj, "try-error")) {
        p.adj <- NA
      }
      category_network$p.adj <- p.adj
    }
    
    rownames(category_network) <- paste(category_network$from,
                                        category_network$to, sep = "-")
  } else {
    category_network <- p.value
  }
  
  return(category_network)
}


#' Fit linear model and compute interaction p-value
#'
#' Helper function to fit a linear model with interaction term and compute
#' the p-value for the interaction coefficient. Used internally by
#' calc_pvalues_network and plot_regressions.
#'
#' @param gene_A_values numeric vector of values for gene A
#' @param gene_B_values numeric vector of values for gene B
#' @param binary_metadata numeric vector (0/1) indicating category membership
#' @return numeric p-value for the interaction term
#' @keywords internal
fit_lm_interaction <- function(gene_A_values, gene_B_values, binary_metadata) {
  design_mat <- cbind(
    1,
    gene_A_values,
    binary_metadata,
    gene_A_values * binary_metadata
  )
  lmfit <- stats::.lm.fit(x = design_mat, y = gene_B_values)
  dof <- length(lmfit$residuals) - length(lmfit$coefficients)
  sigma_squared <- sum((lmfit$residuals)^2) / dof
  XtX_inv <- solve(t(design_mat) %*% design_mat)
  var_covar_matrix <- sigma_squared * XtX_inv
  se_coefficients <- sqrt(diag(var_covar_matrix))
  t_stat <- (lmfit$coefficients[4]) / se_coefficients[4]
  p_interaction <- 2 * (1 - stats::pt(abs(t_stat), df = dof))
  return(p_interaction)
}


#' Get a table of all the significant interactions across categories
#'
#' @param deggs_object an object of class `deggs` generated by
#' `get_diffNetworks`
#' @param assayDataName name of the assayData of interest. If an unnamed list of 
#' data was given to `get_diffNetworks`, `assayDataName` here will be the 
#' number corresponding to the position of the data in the `assayDataList`
#' provided before (i.e. if transcriptomic data was second in the list, a list 
#' of all its differential interactions can be obtained with assayDataName = 2, 
#' if only one data table was provided assayDataName must be 1). Default 1.
#' @param sig_threshold threshold for significance. Default 0.05.
#' @return a `data.frame` listing all the significant differential interactions
#' found across categories for that particular omic data.
#' This list can also be used to substitute or integrate feature selection in 
#' machine learning models for the prediction of the categories (see vignette).
#' @examples
#' data("synthetic_metadata")
#' data("synthetic_rnaseqData")
#' deggs_object <- get_diffNetworks(assayData = synthetic_rnaseqData,
#'                                  metadata = synthetic_metadata,
#'                                  category_variable = "response",
#'                                  verbose = FALSE,
#'                                  show_progressBar = FALSE,
#'                                  cores = 2)
#' get_sig_deggs(deggs_object, sig_threshold = 0.05) 
#' @export
get_sig_deggs <- function(deggs_object, 
                          assayDataName = 1, 
                          sig_threshold = 0.05) {
  
  if (!inherits(deggs_object, "deggs")) (
    stop("Input must be a 'deggs' object")
  )
  
  sig_var <- ifelse(deggs_object[["padj_method"]] == "none", "p.value", "p.adj")
  
  # Extract diffNetworks (=those that are lists) and exlude sig_pvalues_count
  is_list <- vapply(deggs_object[["diffNetworks"]][[assayDataName]],
                    is.list, logical(1))
  diffNetworks_list <- deggs_object[["diffNetworks"]][[assayDataName]][is_list]
  
  if (length(diffNetworks_list) == 0) {
    warning("No valid network found")
    return(data.frame())
  }
  
  # Extract all significant pairs across all categories
  sig_edges_list <- lapply(diffNetworks_list, function(subnetwork) {
    # Check if subnetwork has any rows
    if (is.null(subnetwork) || nrow(subnetwork) == 0) {
      return(data.frame())
    }
    
    # Filter significant edges
    sig_rows <- subnetwork[, sig_var] < sig_threshold
    if (any(sig_rows, na.rm = TRUE)) {
      return(subnetwork[sig_rows, , drop = FALSE])
    } else {
      return(data.frame())
    }
  })
  
  # Remove NULL entries and combine results
  sig_edges_list <- sig_edges_list[!vapply(sig_edges_list, is.null, logical(1))]
  
  if (length(sig_edges_list) == 0) {
    warning("No significant interactions found with ", 
            sig_var, " < ", sig_threshold)
    return(data.frame())
  }
  
  # Combine all significant edges
  sig_edges <- do.call(rbind, sig_edges_list)
  
  # singular column not necessary if mixedModel = FALSE
  if(!deggs_object[["mixedModel"]]) sig_edges$singular <- NULL
  return(sig_edges)
}


#' Get a table of all significant interactions across categories
#'
#' @param deggs_object an object of class `deggs` generated by
#' `get_diffNetworks`
#' @param sig_threshold threshold for significance. Default 0.05.
#' @return a list of multilayer networks (as edge tables), one per 
#' category. 
#' @examples 
#' data("synthetic_metadata")
#' data("synthetic_rnaseqData")
#' data("synthetic_proteomicData")
#' data("synthetic_OlinkData")
#' assayData_list <- list("RNAseq" = synthetic_rnaseqData,
#'                        "Proteomics" = synthetic_proteomicData,
#'                        "Olink" = synthetic_OlinkData)
#' deggs_object <- get_diffNetworks(assayData = assayData_list,
#'                                  metadata = synthetic_metadata,
#'                                  category_variable = "response",
#'                                  verbose = FALSE,
#'                                  show_progressBar = FALSE,
#'                                  cores = 2)
#' get_multiOmics_diffNetworks(deggs_object, sig_threshold = 0.05)
#' @export
get_multiOmics_diffNetworks <- function(deggs_object,
                                        sig_threshold = 0.05) {
  
  if (!inherits(deggs_object, "deggs")) (
    stop("Input must be a 'deggs' object")
  )
  
  if (length(deggs_object[["assayData"]]) == 1) (
    stop("Only one omic dataset was provided. This function is used only in 
         multi-omic scenarios. Use get_sig_deggs() instead.")
  )
  
  categories <- levels(deggs_object[["metadata"]])
  omicDatasets <- names(deggs_object[["assayData"]])
  sig_var <- ifelse(deggs_object[["padj_method"]] == "none", "p.value", "p.adj")
  
  multiLayer_networks <- lapply(categories, function(category) {
    # Collect networks for this category from all omic datasets
    category_networks <- lapply(omicDatasets, function(omicDataset) {
      # Check if the omicDataset entry exists and is a list
      if (!is.list(deggs_object[["diffNetworks"]][[omicDataset]])) {
        return(data.frame())
      }
      
      # Extract the network for this specific category and omicDataset
      network <- deggs_object[["diffNetworks"]][[omicDataset]][[category]]
      
      # Check if the network is a valid data frame
      if (!is.data.frame(network)) {
        return(data.frame())
      }
      
      # Add a source column to identify the omicDataset and filter by sig
      network$layer <- omicDataset
      network <- network[which(network[, sig_var] < sig_threshold), ]
      
      return(network)
    })
    
    # Remove any NULL entries 
    category_networks <- category_networks[!vapply(category_networks, is.null,
                                                   logical(1))]
    
    # Combine all networks for this category
    if (length(category_networks) > 0) {
      merged_network <- do.call(rbind, category_networks)
      merged_network$layer <- as.factor(merged_network$layer)
      # singular column not necessary if mixedModel = FALSE
      if(!deggs_object[["mixedModel"]]) merged_network$singular <- NULL
      return(merged_network)
    } else {
      message("No valid category_networks found for category: ", category, ".")
      return(data.frame())
    }
  })
  names(multiLayer_networks) <- categories
  return(multiLayer_networks)
}

