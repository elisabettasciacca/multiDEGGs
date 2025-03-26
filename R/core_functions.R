#' Generate multi-omic differential networks
#'
#' Generate a multi-layer differential network with interaction p values
#'
#' @param assayData a matrix or data.frame (or list of matrices or data.frames
#' for multi-omic analysis) containing normalised assay data.
#' Sample IDs must be in columns and probe IDs (genes, proteins...) in rows.
#' For multi omic analysis, it is highly recommended to use a named list of data. 
#' If unnamed, sequential names (assayData1, assayData2, etc.) will be
#' assigned to identify each matrix or data.frame. 
#' @param metadata a named vector, matrix, or data.frame containing sample 
#' annotations or categories. If matrix or data.frame, each row should
#' correspond to a sample, with columns representing different sample 
#' characteristics (e.g., treatment group, condition, time point). The colname 
#' of the sample characteristic to be used for differential analysis must be 
#' specified in `category_variable`. Rownames must match the sample IDs used in
#' assayData. 
#' If named vector, each element must correspond to a sample characteristic to 
#' be used for differential andlysis, and names must match sample IDs used in
#' the colnames of `assayData`.
#' Continuous variables are not allowed.  
#' @param category_variable when metadata is a matrix or data.frame this is the 
#' column name of `metadata` that contains the sample annotations to be used for
#' differential analysis
#' @param regression_method whether to use robust linear modelling to calculate
#' link p values. Options are 'rlm' (default) or 'lm'.
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
#' `seq(0.65, 0.98, by = 0.05)`, which generates a sequence of percentiles
#' starting at 0.55, meaning that targets (genes/proteins...) whose expression 
#' value is under the 65th percentile of the whole matrix will be excluded 
#' automatically from the network. Selecting lower starting percentiles allows 
#' more genes (and consequently more links) to be included, but it also 
#' increases computational time.
#' @param use_qvalues whether to use Storey's q values for multiple test
#' adjustment. If FALSE (default), unadjusted p values will be used and shown
#' in the output.
#' @param cores number of cores to use for parallelisation.
#' @return a `deggs` object containing differential networks incorporating
#' p values or adjusted p values for each link.
#' @export
get_diffNetworks <- function(assayData,
                             metadata,
                             category_variable = NULL,
                             regression_method = 'lm',
                             category_subset = NULL,
                             network = NULL,
                             percentile_vector = seq(0.35, 0.98, by = 0.05),
                             use_qvalues = FALSE,
                             show_progressBar = TRUE,
                             verbose = TRUE, 
                             cores = parallel::detectCores() / 2) {
  
  
  if (!(is.list(assayData) || 
        is.matrix(assayData) || is.data.frame(assayData))) (
          stop("assayDataList is not a list or a matrix or data.frame")
        )
  
  if (is.matrix(assayData) || is.data.frame(assayData)) (
    assayDataList <- list(assayData)
  ) else (
    assayDataList <- assayData
  )
  
  if (!is.null(category_variable)) {
    if ((is.matrix(metadata) || is.data.frame(metadata)) && 
        !category_variable %in% colnames(metadata)) {
      stop("category_variable must be %in% colnames(metadata). If metadata is a 
         named vector category_variable must be NULL.")
    }
  } else if (is.matrix(metadata) || is.data.frame(metadata)) {
    stop("metadata is a matrix or dataframe but category_variable is NULL.
       Please specify category_variable to select the correct column 
       from metadata.")
  }
  
  if (is.null(names(assayDataList))) (
    names(assayDataList) <- paste0("assayData", seq_along(assayDataList))
  )
  
  metadata <- tidy_metadata(
    category_subset = category_subset, metadata = metadata,
    category_variable = category_variable, verbose = verbose
  )  # metadata will be a named vector from now on
  
  diffNetworks_list <- lapply(seq_along(assayDataList), function(i){
    assayDataName <- names(assayDataList)[i]
    assayData <- assayDataList[[assayDataName]]
    
    diffNetworks <- 
      tryCatch({
        get_diffNetworks_singleOmic(assayData,
                                    assayDataName, 
                                    metadata,
                                    regression_method,
                                    network,
                                    percentile_vector,
                                    use_qvalues,
                                    show_progressBar,
                                    verbose, 
                                    cores)
      }, error = function(e) { cat(conditionMessage(e)) })
  })
  names(diffNetworks_list) <- names(assayDataList)
  
  degg <- list(
    diffNetworks = diffNetworks_list,
    assayData = assayDataList,
    metadata = metadata,
    regression_method = regression_method,
    category_subset = category_subset,
    use_qvalues = use_qvalues
  )
  class(degg) <- "deggs"
  return(degg)
}

#' Generate differential networks for single onimc analysis 
#' 
#' @param assayData a matrix or data.frame with normalised assay data.
#' Sample IDs must be in columns and probe IDs (genes, proteins...) in rows.
#' @param assayDataName name of the assayData, to identify which omic is. 
#' @param metadata a named vector, each element must correspond to a sample 
#' characteristic to be used for differential analysis, and names must match 
#' sample IDs used in the colnames of `assayData`.
#' Continuous variables are not allowed.  
#' @param regression_method whether to use robust linear modelling to calculate
#' link p values. Options are 'rlm' (default) or 'lm'.
#' @param network network of biological interactions provided by the user. The
#' network must be provided in the form of a table of class data.frame with only 
#' two columns named "from" and "to".
#' If NULL (default) a network of 10,537 molecular interactions obtained from
#' KEGG, mirTARbase, miRecords and transmiR will be used.
#' This has been obtained via the `exportgraph` function of the MITHrIL tool
#' (Alaimo et al., 2016).
#' @param use_qvalues whether to use Storey's q values for multiple test
#' adjustment. If FALSE (default), unadjusted p values will be used and shown
#' in the output.
#' @param cores number of cores to use for parallelisation.
#' @importFrom rlang .data
#' @return a list of differential networks, one per category
get_diffNetworks_singleOmic <- function(assayData,
                                        assayDataName,
                                        metadata,
                                        regression_method = 'rlm',
                                        network = NULL,
                                        percentile_vector = seq(0.55, 0.98,
                                                                by = 0.05), 
                                        use_qvalues = FALSE,
                                        show_progressBar = TRUE,
                                        verbose = TRUE, 
                                        cores = parallel::detectCores() / 2) {
  
  if (verbose) (
    message(paste0("\nGenerating ", assayDataName, " network layer..."))
  )
  sig_var <- ifelse(use_qvalues, "q.value", "p.value")
  
  if (!(is.matrix(assayData) || is.data.frame(assayData))) (
    stop(paste0(assayDataName, " is not a matrix or data.frame"))
  )
  
  if(!any(colnames(assayData) %in% names(metadata))) (
    stop(paste0("Sample IDs in ", assayDataName, " don't match the IDs in 
                metadata."))
  )
  
  # align metadata and assayData
  if (length(intersect(colnames(assayData), names(metadata))) == 0) {
    stop(paste0("Sample IDs in metadata and ", assayDataName, " don't match."))
  }
  
  # remove assayData samples that don't exist in metadata (because they
  # were missing or because they were filtered out due to the category_subset)
  assayData <- assayData[, which(colnames(assayData) %in% names(metadata)),
                         drop = FALSE]
  
  # check and remove metadata sample IDs that don't exist in assayData
  missing_metadataSamples <- names(metadata)[which(!(names(metadata) %in% 
                                                       colnames(assayData)))]
  
  if (length(missing_metadataSamples) != 0) {
    missing_samples_formatted <- paste(missing_metadataSamples, collapse = ", ")
    message(paste0("\nThe following samples IDs are missing in ", assayDataName, 
                   ":\n", missing_samples_formatted))
    metadata <- metadata[!(names(metadata) %in% missing_metadataSamples)]
    
    if(length(unique(metadata)) == 1) (
      stop(paste0("All sample IDs in ", assayDataName, " belong to one 
                  category. No differential analysis is possible."))
    )
  }
  
  if (!is.null(network)) (
    # user provided network
    network_to_use <- network
  ) else (
      network_to_use <- metapathway_gene_symbols
  )
  
  edges <- network_to_use %>%
    dplyr::filter(.data$from %in% rownames(assayData)) %>%
    dplyr::filter(.data$to %in% rownames(assayData))  
  
  if (nrow(edges) == 0) (
    stop(paste0("Rownames of ", assayDataName, " don't match any link of the 
                biological network. Check if names need conversion to gene 
                symbols or entrez IDs. Alternatively, provide more data."))
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
    category_vec <- apply(category_vec, 1, median, na.rm = TRUE)
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
      "assayData", "calc_pvalues_network2", "metadata",
      "sig_edges_count", "sig_var"
    ), envir = environment())
    
    parallel::clusterEvalQ(cl, {
      library("dplyr")
    })
    if (show_progressBar) (
      pvalues_list <- pbapply::pblapply(
        cl = cl, percentile_vector,
        function(percentile) {
          calc_pvalues_percentile(
            assayData = assayData,
            sig_var = sig_var,
            metadata = metadata,
            percentile = percentile,
            category_median_list = category_median_list,
            contrasts = contrasts,
            regression_method = regression_method,
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
            sig_var = sig_var,
            metadata = metadata,
            percentile = percentile,
            category_median_list = category_median_list,
            contrasts = contrasts,
            regression_method = regression_method,
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
            sig_var = sig_var,
            metadata = metadata,
            percentile = percentile,
            category_median_list = category_median_list,
            contrasts = contrasts,
            regression_method = regression_method,
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
            sig_var = sig_var,
            metadata = metadata,
            percentile = percentile,
            category_median_list = category_median_list,
            contrasts = contrasts,
            regression_method = regression_method,
            edges = edges,
            categories_length = length(categories),
            sig_edges_count = sig_edges_count
          )
        )
      )
    )
    
  }
  names(pvalues_list) <- percentile_vector
  
  # Extract the network filtered at the optimal threshold percentile
  # (the one with the highest number of statistically significant interactions)
  sig_pvalues <- sapply(pvalues_list, function(networks) {
    if (!inherits(networks, "try-error") && !is.null(networks)) {
      return(networks$sig_pvalues_count)
    }
  }) %>% unlist()
  
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
  return(final_networks)
}


#' Tidying up of metadata. Samples belonging to undesidered categories
#' (if specified) will be removed as well as categories with less than five 
#' samples, and NAs.
#'
#' @param category_subset optional character vector indicating which categories
#' are used for comparison. If not specified, all categories in
#' `category_variable` will be used.
#' @param metadata a data.frame, matrix, or named vector of annotations.
#' @param category_variable column name in `metadata` (if data.frame or matrix)
#' or NULL if `metadata` is already a named vector containing category 
#' information.
#' @param verbose logical; whether to print messages about processing steps.
#' @return a tidy named factor vector of sample annotations.
tidy_metadata <- function(category_subset = NULL,
                          metadata,
                          category_variable = NULL, 
                          verbose = FALSE) {
  
  # Convert metadata to a named vector based on input type
  if (is.atomic(metadata) || is.factor(metadata)) {
    if (is.null(names(metadata))) {
      stop("Named vector required: metadata vector must have names.")
    }
    # Already a named vector
    category_vector <- metadata
  } else if (is.matrix(metadata) || is.data.frame(metadata)) {
    if (is.null(category_variable)) {
      stop("category_variable must be specified when metadata is a matrix or data.frame.")
    }
    if (!category_variable %in% colnames(metadata)) {
      stop("category_variable '", category_variable, "' not found in metadata.")
    }
    category_vector <- metadata[, category_variable]
    names(category_vector) <- rownames(metadata)
  } else {
    stop("metadata must be a named vector/factor, matrix, or data.frame.")
  }
  
  # Subset seleted categories (optional)
  if (!is.null(category_subset)) {
    if (!all(category_subset %in% category_vector)) (
      stop("category_subset values must be contained in metadata")
    )
    keep_samples <- names(category_vector)[category_vector %in% category_subset]
    if (length(keep_samples) == 0) {
      stop("No samples remain after filtering for specified category_subset.")
    }
    category_vector <- category_vector[keep_samples]
  }
  
  # Convert to factor
  if (!is.factor(category_vector)) {
    category_vector <- as.factor(category_vector)
    if (verbose) {
      message("category values were converted to factor.")
    }
  }
  
  # Remove categories with less than five observations
  tbl <- table(category_vector)
  small_groups <- names(tbl)[tbl < 5]
  
  if (length(small_groups) > 0) {
    for (group in small_groups) {
      message("The ", group, 
              " category did not contain enough samples (less than five observations). ",
              "This category will be removed.")
    }
    keep_samples <- names(category_vector)[category_vector %in% names(tbl)[tbl >= 5]]
    category_vector <- category_vector[keep_samples]
    category_vector <- droplevels(category_vector)
    
    if (nlevels(category_vector) < 2) {
      stop("At least two categories with at least 5 observations are required.")
    }
  }
  
  # Remove NAs
  NAs_number <- sum(is.na(category_vector))
  if (NAs_number > 0) {
    category_vector <- na.omit(category_vector)
    if (verbose) {
      message("category values contained NAs. ", NAs_number, " samples were removed.")
    }
  }
  return(category_vector)
}


#' Compute interaction p values for a single percentile value
#'
#' @inheritParams get_diffNetworks
#' @param categories_length integer number indicating the number of categories
#' @param category_median_list list of category data.frames
#' @param sig_var Inherited from `get_diffNetworks`. It can be
#' `q.value` or `p.value` depending on how `use_qvalues` was set in the
#' `get_diffNetworks` function (default `FALSE`).
#' @param percentile a float number indicating the percentile to use.
#' @param contrasts data.frame containing the categories contrasts in rows
#' @param regression_method whether to use robust linear modelling to obtain
#' p value of the interactions. Options are 'rlm' (default) or 'lm'
#' @param edges network of biological interactions in the form of a table of
#' class data.frame with two columns: "from" and "to".
#' @param sig_edges_count number of significant edges (p < 0.05)
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @return The list of float numbers of the significant pvalues
#' for a specific percentile
calc_pvalues_percentile <- function(assayData,
                                    metadata,
                                    categories_length,
                                    category_median_list,
                                    sig_var,
                                    percentile,
                                    contrasts,
                                    regression_method = "rlm",
                                    edges,
                                    sig_edges_count) {
  # 1st filtering step: remove low expressed genes in each category
  # (genes under the percentile threshold)
  cut_off <- stats::quantile(as.matrix(assayData), prob = percentile, na.rm = T)
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
      categoryEdges <- edges %>%
        dplyr::filter(.data$from %in% names(category_vec)) %>%
        dplyr::filter(.data$to %in% names(category_vec))
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
  tot_edges <- unlist(lapply(categories_network_list, function(category_network) {
    if (!is.character(category_network)) (
      nrow(category_network)
    )
  })) %>% sum()
  
  # calculate interaction p values
  pvalues_list <- lapply(categories_network_list, function(category_network) {
    return(calc_pvalues_network2(
      category_network = category_network,
      assayData = assayData,
      metadata = metadata,
      regression_method = regression_method,
      categories_length = categories_length,
      sig_var = sig_var
    ))
  })
  # count significant p values
  if (tot_edges > sig_edges_count) {
    # num tot edges greater than previous sig edges count
    p_values_sig_count <- unlist(lapply(pvalues_list, function(category_network) {
      if (!is.null(dim(category_network))) {
        num_sig_links <- sum(category_network[, sig_var] < 0.05, na.rm = TRUE)
        # return the number of links with significant pvalues or qvalues
        if (!is.na(num_sig_links)) {
          return(num_sig_links)
        } else {return(0)}
      } else {return(0)}
    }))
    num_sig_pvalues <- sum(p_values_sig_count) # these are pvalues or qvalues
    
    # The <<- operator here is used only for efficiency reasons in a controlled 
    # fashion.
    # <<- represents the only way to implement early stopping of the loop 
    # when the total number of edges left in the network for a given 
    # percentile is smaller that the maximum number of significant edges found 
    # in previous percentile's networks. It wouldn't make sense to keep looking 
    # for a bigger number of significant edges in subsequent percentiles. 
    # The update of sig_edges_count will allow skipping unnecessary lm or rlm
    # iterations. 
    if (num_sig_pvalues > sig_edges_count) {
      sig_edges_count <<- num_sig_pvalues
    }
    pvalues_list <- append(pvalues_list, num_sig_pvalues)
    names(pvalues_list)[length(pvalues_list)] <- "sig_pvalues_count"
    return(pvalues_list)
  } else {
    # previous sig edges count greater than num tot edges, skip all calculations
    return(NULL)
  }
}

#' Calculate the pvalues for specific category network samples
#'
#' @inheritParams calc_pvalues_percentile
#' @param category_network network table for a specific category
#' @importFrom methods is
#' @return a list of p values
calc_pvalues_network <- function(assayData,
                                 metadata,
                                 sig_var,
                                 categories_length,
                                 regression_method = 'rlm',
                                 category_network) {
  
  if (!is.character(category_network)) {
    if (nrow(category_network) > 0) {
      # prepare data
      df_list <- mapply(function(gene_B, gene_A) {
        # using both t() and as.vector to be compatible with both matrices and dfs
        df <- data.frame(as.vector(t(assayData[gene_A, ])),
                         as.vector(t(assayData[gene_B, ])),
                         metadata,
                         check.names = FALSE)
        colnames(df)[1:2] <- c(gene_A, gene_B)
        return(df)
      }, gene_A = category_network$from, gene_B = category_network$to,
      SIMPLIFY = FALSE)
      
      # calculate regressions with interaction term
      p_values <- lapply(df_list, function(df) {
        if (categories_length == 2) {
          if (regression_method == "lm") {
            # gene_B ~ gene_A * category
            lmfit <- stats::lm(df[, 2] ~ df[, 1] * df[, 3])
            p_interaction <- stats::coef(summary(lmfit))[4, 4]
          }
          if (regression_method == "rlm") {
            # gene_B ~ gene_A * category
            robustfit <- suppressWarnings(MASS::rlm(df[, 2] ~ df[, 1] * df[, 3]))
            p_interaction <- try(sfsmisc::f.robftest(robustfit, var = 3)$p.value,
                                 silent = TRUE)
            if (class(p_interaction) == "try-error") (
              p_interaction <- NA
            )
          }
          output <- data.frame(
            from = colnames(df)[1], to = colnames(df)[2],
            p.value = p_interaction
          )
        }
        if (categories_length >= 3) {
          # one-way ANOVA
          # gene_B ~ gene_A * category
          res_aov <- stats::aov(df[, 2] ~ df[, 1] * df[, 3], data = df)
          p_interaction <- summary(res_aov)[[1]][["Pr(>F)"]][3]
          output <- data.frame(
            from = colnames(df)[1], to = colnames(df)[2],
            p.value = p_interaction
          )
        }
        return(output)
      })
    } else {
      p_values <- "No specific links for this category."
    }
  } else {
    p_values <- "No specific links for this category."
  }
  
  if (is.list(p_values)) {
    # make a data.frame with all values
    p_values <- as.data.frame(do.call("rbind", p_values))
    p_values$from <- as.character(p_values$from)
    p_values$to <- as.character(p_values$to)
    
    if (sig_var == "q.value") {
      # adding Storey's q values
      q.values <- try(qvalue::qvalue(p_values[, "p.value"])$qvalues,
                      silent = TRUE)
      if (is(q.values, "try-error")) {
        if (nrow(p_values) > 1) (
          q.values <- qvalue::qvalue(p = p_values[, "p.value"], pi0 = 1)$qvalues
        ) else (
          q.values <- NA
        )
      }
      p_values$q.value <- q.values
    }
  }
  if (!is.character(p_values)) {
    rownames(p_values) <- paste(p_values$from, p_values$to, sep = "-")
  }
  return(p_values)
}


#' Calculate the pvalues for specific category network samples
#'
#' @inheritParams calc_pvalues_percentile
#' @param category_network network table for a specific category
#' @importFrom methods is
#' @return a list of p values
calc_pvalues_network2 <- function(assayData,
                                 metadata,
                                 sig_var,
                                 categories_length,
                                 regression_method = 'lm',
                                 category_network) {
  
  if (!is.character(category_network)) {
    if (nrow(category_network) > 0) {
      p.value <- sapply(1:nrow(category_network), function(i){
        gene_A <- category_network[i, 1]
        gene_B <- category_network[i, 2]
        
        if (categories_length == 2) {
          if (regression_method == "lm") {
            binary.metadata <- as.numeric(metadata) - 1
            design.mat <- cbind(
              1,
              as.numeric(assayData[gene_A, ]),
              binary.metadata,
              as.numeric(assayData[gene_A, ]) * binary.metadata
              )
            lmfit <- stats::.lm.fit(x = design.mat,
                                    y = as.numeric(assayData[gene_B, ]))
            dof <- length(lmfit$residuals) - length(lmfit$coefficients)
            # Mean Squared Error
            sigma_squared <- sum((lmfit$residuals)^2) / dof
            # variance-covariance of coefficients
            XtX_inv <- solve(t(design.mat) %*% design.mat)
            var_covar_matrix <- sigma_squared * XtX_inv
            se_coefficients <- sqrt(diag(var_covar_matrix))
            # t-statics for the interaction term (4th coeff) 
            t_stat <- (lmfit$coefficients[4]) / se_coefficients[4]
            p_interaction <- 2 * (1 - pt(abs(t_stat), df = dof))
          }
          if (regression_method == "rlm") {
            # gene_B ~ gene_A * category
            robustfit <- suppressWarnings(
              MASS::rlm(as.numeric(assayData[gene_B, ]) ~
                          as.numeric(assayData[gene_A, ]) * metadata)
              )
            p_interaction <- try(
              sfsmisc::f.robftest(robustfit, var = 3)$p.value, silent = TRUE
              )
            if (class(p_interaction) == "try-error") (
              p_interaction <- NA
            )
          }
        }
        if (categories_length >= 3) {
          # one-way ANOVA
          # gene_B ~ gene_A * category
          res_aov <- stats::aov(as.numeric(assayData[gene_B, ]) ~
                                  as.numeric(assayData[gene_A, ]) * metadata)
          p_interaction <- summary(res_aov)[[1]][["Pr(>F)"]][3]
        }
        return(p_interaction)
      })
    } else {
      p.value <- "No specific links for this category."
    }
  } else {
    p.value <- "No specific links for this category."
  }
  
  if (is.numeric(p.value)) {
    category_network <- cbind(category_network, p.value)
    
    if (sig_var == "q.value") {
      # adding Storey's q values
      q.values <- try(qvalue::qvalue(category_network[, "p.value"])$qvalues,
                      silent = TRUE)
      if (is(q.values, "try-error")) {
        if (nrow(category_network) > 1) (
          q.values <- qvalue::qvalue(p = category_network[, "p.value"], pi0 = 1)$qvalues
        ) else (
          q.values <- NA
        )
      }
      category_network$q.value <- q.values
    }
    rownames(category_network) <- paste(category_network$from,
                                        category_network$to, sep = "-")
  } else {
    category_network <- p.value
  }
  return(category_network)
}

#' Get a table of all the significant interactions across categories
#'
#' @param deggs_object an object of class `deggs` generated by
#' `get_diffNetworks`
#' @param assayDataName name of the assayData of interest. If an unnamed list of 
#' data was given to `get_diffNetworks`, the assayDataName here will be the 
#' number indicating the position of the data in the assayDataList provided
#' before (i.e. if transcriptomic data was second in the list, a list of 
#' all its differential interactions can be obtained with assayDataName = 2, 
#' if only one data table was provided assayDataName must be 1). Default 1.
#' @param sig_threshold threshold for significance. Default 0.05.
#' @return a `data.frame` listing all the significant differential interactions
#' found across categories for that particular omic data.
#' This list can also be used to substitute or integrate feature selection in 
#' machine learning models for the prediction of the categories.
#' @export
get_sig_deggs <- function(deggs_object, 
                          assayDataName = 1, 
                          sig_threshold = 0.05) {
  
  if (!inherits(deggs_object, "deggs")) (
    stop("Input must be a 'deggs' object")
  )
  
  sig_var <- ifelse(deggs_object[["use_qvalues"]], "q.value", "p.value")
  
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
  sig_edges_list <- sig_edges_list[!sapply(sig_edges_list, is.null)]
  
  if (length(sig_edges_list) == 0) {
    warning("No significant interactions found with ",
            sig_var, " < ", sig_threshold)
    return(data.frame())
  }
  
  # Combine all significant edges
  sig_edges <- do.call(rbind, sig_edges_list)
  return(sig_edges)
}


#' Get a table of all significant interactions across categories
#'
#' @param deggs_object an object of class `deggs` generated by
#' `get_diffNetworks`
#' @param sig_threshold threshold for significance. Default 0.05.
#' @return a list of multilayer networks (in form of edge table), one per 
#' category in the data. 
#' @export
get_multiOmics_diffNetworks <- function(deggs_object,
                                        sig_threshold = 0.05) {
  
  if (!inherits(deggs_object, "deggs")) (
    stop("Input must be a 'deggs' object")
  )
  
  if (length(deggs_object[["assayData"]]) == 1) (
    stop("Only one omic dataset was provided. This function is used only in 
         multi-omic scenarios.")
  )
  
  categories <- levels(deggs_object[["metadata"]])
  omicDatasets <- names(deggs_object[["assayData"]])
  sig_var <- ifelse(deggs_object[["use_qvalues"]], "q.value", "p.value")
  
  multiLayer_networks <- lapply(categories, function(category) {
    # Collect networks for this category from all omic datasets
    category_networks <- lapply(omicDatasets, function(omicDataset) {
      # Check if the omicDataset entry exists and is a list
      if (!is.list(deggs_object[["diffNetworks"]][[omicDataset]])) {
        warning(paste("Dataset", omicDataset, "is not a list. Skipping."))
        return(data.frame())
      }
      
      # Extract the network for this specific category and omicDataset
      network <- deggs_object[["diffNetworks"]][[omicDataset]][[category]]
      
      # Check if the network is a valid data frame
      if (!is.data.frame(network)) {
        warning(paste("Category", category, "in dataset", omicDataset, 
                      "is not a data frame. Skipping."))
        return(data.frame())
      }
      
      # Add a source column to identify the omicDataset
      network$layer <- omicDataset

      return(network)
    })
    
    # Remove any NULL entries 
    category_networks <- category_networks[!sapply(category_networks, is.null)]
    
    # Combine all networks for this category
    if (length(category_networks) > 0) {
      merged_network <- do.call(rbind, category_networks)
      merged_network$layer <- as.factor(merged_network$layer)
      rownames(merged_network) <- paste0(rownames(merged_network), "_", 
                                         merged_network$layer)
      merged_network <- merged_network[which(merged_network[, sig_var] < 
                                               sig_threshold), ]
      return(merged_network)
    } else {
      warning(paste("No valid category_networks found for category:", category))
      return(data.frame())
    }
  })
  names(multiLayer_networks) <- categories
  return(multiLayer_networks)
}

