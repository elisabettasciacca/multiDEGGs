#' multiDEGGs_filter
#' 
#' Function to be passed to the `modifyX` parameter of 
#' \link[nestedcv]{nestcv.train} or \link[nestedcv]{nestcv.glmnet} to allow 
#' nested feature selection and augmentation via differential network analysis 
#' with multiDEGGs.
#' 
#' @param y Numeric vector or factor. Response variable (outcome), i.e. 
#'  the 'metadata' named vector, as passed by \link[nestedcv]{nestcv.train} or 
#'  \link[nestedcv]{nestcv.glmnet}.
#' @param x Predictor variables, i.e. the assayData matrix with genes in columns
#'  and IDs in rows, as passed by \link[nestedcv]{nestcv.train} or
#'   \link[nestedcv]{nestcv.glmnet}.
#' @param keep_single_genes Logical, default FALSE. If TRUE, the function will 
#' return unique individual genes along with significant pairs. 
#' @param nfilter Integer. Maximum total number of predictors to return. 
#' When `keep_single_genes = TRUE`, this parameter limits the combined count of 
#' unique and paired predictors
#' (i.e., length(keep_DEGGs) + nrow(pairs) <= nfilter). 
#' Predictors are included from most to least significant: for each pair, both 
#' the pair itself and the new unique variables are included until the `nfilter`
#' threshold is reached. 
#' When `keep_single_genes = FALSE`, `nfilter` only limits the number of pairs 
#' returned. 
#' If `NULL`, no filtering is applied and the total number of predictors 
#' will depend on how many significantly different interactions are detected by
#' multiDEGGs in that fold.
#'  
#' @returns a list containing two types of predictors:
#'  - single predictors (stored in the 'keep_DEGGs' vector) 
#'  - paired predictors (stored in the 'pairs' dataframe)
#' Note that `nfilter` limits the maximum number of engineered features returned,
#' however this number might be lower and will depend on how many significantly 
#' different interactions will be found in each fold by multiDEGGs.
#' \bold{If no significantly different interactions are found
#' the function will print a '0' and switch to t-test for that fold.}
#' 
#' @examples
#' library(nestedcv)
#' data("synthetic_metadata")
#' data("synthetic_rnaseqData")
#' 
#' # fit a regularized linear model
#' # Note that nfilter, n_outer_folds, n_inner_folds are set low to keep the
#' # example lightweight. Adjust these values as needed for your use case.
#' \dontrun{
#' fit.glmnet <- nestcv.glmnet(
#'   y = as.numeric(synthetic_metadata$response),
#'   x =  t(synthetic_rnaseqData),
#'   modifyX = "multiDEGGs_filter",
#'   modifyX_options = list(keep_single_genes = FALSE,
#'                          nfilter = 5),
#'   modifyX_useY = TRUE,
#'   n_outer_folds = 4,
#'   n_inner_folds = 4)
#' 
#' summary(fit.glmnet)
#' }
#' 
#' # fit a random forest model:
#' # note that nfilter, n_outer_folds, n_inner_folds are set low to keep the
#' # example lightweight. Adjust these values as needed for your use case.
#' fit.rf <- nestcv.train(
#'   y = synthetic_metadata$response,
#'   x = t(synthetic_rnaseqData),
#'   method = "rf",
#'   modifyX = "multiDEGGs_filter",
#'   modifyX_options = list(keep_single_genes = FALSE,
#'                          nfilter = 5),
#'   modifyX_useY = TRUE,
#'   n_outer_folds = 2,
#'   n_inner_folds = 2
#' )
#' 
#' fit.rf$summary
#' @export
multiDEGGs_filter <- function(y,
                              x, 
                              keep_single_genes = FALSE,
                              nfilter = NULL) {
  counts <- t(x)
  names(y) <- colnames(counts)
  deggs_object <- get_diffNetworks(assayData = counts,
                                   metadata = y,
                                   regression_method = "lm",
                                   percentile_vector = seq(0.25, 0.98, by = 0.05),
                                   padj_method = "q.value", 
                                   show_progressBar = FALSE,
                                   verbose = FALSE,
                                   cores = 1)
  
  pairs <- get_sig_deggs(deggs_object, 1, 0.05)
  
  if(keep_single_genes) {
    if (!is.null(nfilter)) {
      # Process pairs in order of significance 
      # and collect elements until nfilter is reached
      keep_DEGGs <- c()
      selected_pairs_idx <- c()
      total_elements <- 0
      
      for (i in 1:nrow(pairs)) {
        if (total_elements >= nfilter) break
        
        # Get genes from current pair
        current_genes <- c(pairs[i, 1], pairs[i, 2])
        new_genes <- setdiff(current_genes, keep_DEGGs)
        
        # Check if we can add this pair (1 element) plus any new genes
        elements_to_add <- 1 + length(new_genes)  # 1 pair + new genes
        
        if (total_elements + elements_to_add <= nfilter) {
          # Add the pair and new genes
          selected_pairs_idx <- c(selected_pairs_idx, i)
          keep_DEGGs <- c(keep_DEGGs, new_genes)
          total_elements <- total_elements + elements_to_add
        } else {
          # Check if we can add only the pair (if both genes already exist)
          if (length(new_genes) == 0 && total_elements + 1 <= nfilter) {
            selected_pairs_idx <- c(selected_pairs_idx, i)
            total_elements <- total_elements + 1
          } else {
            # Check if we can add only some new genes (without the pair)
            remaining_slots <- nfilter - total_elements
            if (remaining_slots > 0 && length(new_genes) > 0) {
              genes_to_add <- new_genes[1:min(remaining_slots,
                                              length(new_genes))]
              keep_DEGGs <- c(keep_DEGGs, genes_to_add)
              total_elements <- total_elements + length(genes_to_add)
            }
            break
          }
        }
      }
      
      # Filter pairs to only selected indices
      pairs <- pairs[selected_pairs_idx, , drop = FALSE]
    } else {
      # No nfilter applied, extract all unique genes from pairs
      keep_DEGGs <- unique(unlist(lapply(1:nrow(pairs), function(i)(
        row_n = c(pairs[i,1], pairs[i,2])
      ))))
    }
  } else {
    keep_DEGGs <- c()
    if (!is.null(nfilter)) {
      # When keep_single_genes is FALSE, nfilter only limits pairs
      if(nrow(pairs) > nfilter){
        pairs <- pairs[1:nfilter, ]
      }
    }
  }

  selected_features <- list(keep = keep_DEGGs, pairs = pairs)

  # If no differential interaction is found, default to t-test filter
  if(nrow(selected_features$pairs) <= 1){
    cat_parallel("0")
    keep_ttest <- nestedcv::ttest_filter(y, x, nfilter = nfilter, 
                                         type = "names")
    selected_features <- list(keep = keep_ttest, pairs = data.frame())
  }
  
  class(selected_features) <- "multiDEGGs_filter"
  return(selected_features)
}


##' Combined multiDEGGs filter
#'
#' This function can be passed to the `modifyX` parameter of 
#' \link[nestedcv]{nestcv.train} or \link[nestedcv]{nestcv.glmnet}
#' to use one of the available statistical filters (t-test, wilcoxon, etc.) in 
#' combination with multiDEGGs. 
#' Single predictors will be selected by the selected statistical filter an 
#' paired predictors will be added by multiDEGGs. 
#' 
#' @param y Numeric vector or factor. Response variable (outcome), i.e. 
#'  the 'metadata' named vector, as passed by \link[nestedcv]{nestcv.train} or 
#'  \link[nestedcv]{nestcv.glmnet}.
#' @param x Predictor variables, i.e. the assayData matrix with genes in columns
#'  and IDs in rows, as passed by \link[nestedcv]{nestcv.train} or
#'   \link[nestedcv]{nestcv.glmnet}.
#' @param filter_method Character string. Statistical filtering method to be
#' used in combination with multiDEGGs for sigle feature selection. 
#' Options are: "ttest", "wilcoxon", "ranger", "glmnet", "pls".
#' @param nfilter Integer. Maximum number of features to select. 
#' @param dynamic_nfilter Logical. If `TRUE` `nfilter` will limit the number of 
#' features selected by the statistical filter and the feature space will be
#' augmented by adding ALL the paired predictors found by multiDEGGs.
#' If `FALSE` `nfilter` will limit the total number of predictors, with 
#' approximately half allocated to pairs and half to single genes.
#' @param keep_single_genes Logical. When `dynamic_nfilter = TRUE`, determines
#'   whether to include single genes selected by multiDEGGs (i.e. the single 
#'   variables included in the differential pairs) in addition to those from
#'   the statistical filter. Default is FALSE.
#' @param ... Additional arguments passed to the filtering functions.
#'
#' @details
#' The function operates in two modes:
#'
#' \strong{Dynamic Filtering (dynamic_nfilter = TRUE):}
#' - Selects \code{nfilter} single genes using the specified statistical method
#' - Finds all significant gene pairs using multiDEGGs
#' - Total predictors = nfilter single genes + number of significant pairs
#'    - If \code{keep_single_genes = TRUE}, also includes single genes obtained 
#'      from pairs found by multiDEGGs
#'
#' \strong{Balanced Selection (dynamic_nfilter = FALSE):}
#' - Allocates approximately half of \code{nfilter} to gene pairs
#' - Remaining slots filled with single genes from the statistical filter
#' - If fewer pairs are found than allocated, compensates by selecting more 
#'   single genes
#' - Ensures consistent total number of predictors across outer folds
#'
#' The statistical filtering methods include:
#' - \code{"ttest"}: Two-sample t-test for differential expression
#' - \code{"wilcoxon"}: Wilcoxon rank-sum test
#' - \code{"ranger"}: Random Forest variable importance
#' - \code{"glmnet"}: Elastic net regularization
#' - \code{"pls"}: Partial Least Squares variable importance
#'
#' @return An object of class "multiDEGGs_filter" containing:
#' \item{keep}{Character vector of selected single gene names}
#' \item{pairs}{Data frame of selected gene pairs with interaction information}
#' @examples
#' library(nestedcv)
#' data("synthetic_metadata")
#' data("synthetic_rnaseqData")
#' 
#' # fit a regularized linear model
#' # note that nfilter, n_outer_folds, n_inner_folds are set low to keep the
#' # example lightweight. Adjust these values as needed for your use case.
#' \dontrun{
#' fit.glmnet <- nestedcv::nestcv.glmnet(
#' y = as.numeric(synthetic_metadata$response),
#' x =  t(synthetic_rnaseqData),
#' modifyX = "multiDEGGs_combined_filter",
#' modifyX_options = list(filter_method = "ttest", 
#'                        nfilter = 5,
#'                        dynamic_nfilter = TRUE, 
#'                        keep_single_genes = FALSE),
#' modifyX_useY = TRUE,
#' n_outer_folds = 4,
#' n_inner_folds = 4)
#' 
#' summary(fit.glmnet)
#' }
#' 
#' # fit a random forest model
#' # NOTE: nfilter, n_outer_folds, n_inner_folds are set low to keep the
#' # example lightweight. Adjust these values as needed for your use case.
#' fit.rf <- nestedcv::nestcv.train(
#'   y = synthetic_metadata$response,
#'   x = t(synthetic_rnaseqData),
#'   method = "rf",
#'   modifyX = "multiDEGGs_combined_filter",
#'   modifyX_options = list(filter_method = "ttest", 
#'                          nfilter = 5,
#'                          dynamic_nfilter = TRUE, 
#'                          keep_single_genes = FALSE),
#'   modifyX_useY = TRUE,
#'   n_outer_folds = 2,
#'   n_inner_folds = 2
#' )
#' 
#' fit.rf$summary
#' 
#' @export
multiDEGGs_combined_filter <- function(y,
                                       x,
                                       filter_method = "ttest", 
                                       nfilter,
                                       dynamic_nfilter = TRUE, 
                                       keep_single_genes = FALSE, 
                                       ...) {
  if (dynamic_nfilter) {
    # The total number of predictors will depend on how many sig interactions 
    # are found by multiDEGGs, so the tot number of predictors will be 
    # nfilter + pairs in each outer fold
    keep_filtered <- switch(filter_method,
                            "ttest" = nestedcv::ttest_filter(y, x, nfilter = nfilter,
                                                   type = "names"),
                            "wilcoxon" = nestedcv::wilcoxon_filter(y, x, nfilter = nfilter,
                                                         type = "names"),
                            "ranger" = nestedcv::ranger_filter(y, x, nfilter = nfilter, 
                                                     type = "names"),
                            "glmnet" = nestedcv::glmnet_filter(y, x, nfilter = nfilter, 
                                                     type = "names"),
                            "pls" = nestedcv::pls_filter(y, x, nfilter = nfilter, 
                                               type = "names"),
                            stop("Invalid filter_method. 
                                 Choose from: ttest, wilcoxon, ranger, 
                                              glmnet, pls"))
    
    multiDEGGs_res <- multiDEGGs_filter(y, x, keep_single_genes = TRUE)

    if (keep_single_genes) {
      # select single genes from both multiDEGGs and chosen filter
      selected_features <- list(keep = c(keep_filtered, 
                                         multiDEGGs_res$keep_DEGGs), 
                                pairs = multiDEGGs_res$pairs)
    } else {
      # single predictors from filter_method, pairs from multiDEGGs
      selected_features <- list(keep = keep_filtered,
                                pairs = multiDEGGs_res$pairs)
    }
  } else {
    # balanced selection, no single genes selected by multiDEGGs (only pairs)
    n_pairs <- floor(nfilter / 2)   # max half for interactions
    n_single <- nfilter - n_pairs   # remaining for single genes
    
    # Get significant pairs from multiDEGGs first (no filter on number)
    multiDEGGs_res <- multiDEGGs_filter(y, x, nfilter = NULL,
                                        keep_single_genes = FALSE)
    
    # Select available pairs (up to n_pairs)
    if (nrow(multiDEGGs_res$pairs) > 0) {
      selected_pairs <- multiDEGGs_res$pairs[1:min(n_pairs, 
                                                   nrow(multiDEGGs_res$pairs)),,
                                             drop = FALSE]
      actual_pairs <- nrow(selected_pairs)
    } else {
      selected_pairs <- data.frame()
      actual_pairs <- 0
    }
    
    # Calculate how many single genes we need from selected filter 
    # (compensate for missing pairs to have an even number of predictors in each
    # outer fold)
    singles_needed <- n_single + (n_pairs - actual_pairs)
    
    # Get single genes from selected filter
    keep_filtered <- switch(filter_method,
                            "ttest" = nestedcv::ttest_filter(y, x,
                                                   nfilter = singles_needed,
                                                   type = "names"),
                            "wilcoxon" = nestedcv::wilcoxon_filter(y, x, 
                                                    nfilter = singles_needed,
                                                    type = "names"),
                            "ranger" = nestedcv::ranger_filter(y, x, 
                                                     nfilter = singles_needed, 
                                                     type = "names"),
                            "glmnet" = nestedcv::glmnet_filter(y, x, 
                                                     nfilter = singles_needed,
                                                     type = "names"),
                            "pls" = nestedcv::pls_filter(y, x, 
                                               nfilter = singles_needed, 
                                               type = "names"),
                            stop("Invalid filter_method. 
                                 Choose from: ttest, wilcoxon, 
                                 ranger, glmnet, pls")
    )
    selected_features <- list(keep = keep_filtered, pairs = selected_pairs)
  }
  class(selected_features) <- "multiDEGGs_filter"
  return(selected_features)
}


#' Predict method for multiDEGGs_filter objects
#'
#' This function generates predictions by creating a dataset with single 
#' and combined predictors based on the filtering results of a
#' multiDEGGs_filter model.
#'
#' @param object A fitted object of class `multiDEGGs_filter` containing
#'  filtering results with:
#'   \describe{
#'     \item{keep}{Character vector of variable names to keep as single 
#'     predictors}
#'     \item{pairs}{Data frame or matrix with two columns specifying pairs 
#'     of variables to combine}
#'   }
#' @param newdata A data frame containing the new data for prediction. 
#' Must contain all variables specified in \code{object$keep} and
#' \code{object$pairs}.
#' @param interaction.type Character string specifying how to combine the 
#' paired predictors.
#'   Options are:
#'   \describe{
#'     \item{"ratio"}{Combine paired predictors by dividing the first variable
#'      by the second (a/b)}
#'     \item{other}{Combine paired predictors by multiplying the variables (a*b)
#'     }
#'   }
#'   Default is "ratio".
#' @param sep Character string used as separator when creating column names for
#'   combined predictors. Default is ":".
#' @param ... Additional arguments passed to the generic function.
#'
#' @return A data frame containing:
#'   \itemize{
#'     \item Single predictors (if any are specified in \code{object$keep})
#'     \item Combined predictors based on variable pairs and interaction type
#'   }
#'
#' @details
#' The function processes the filtering results in two steps:
#' \enumerate{
#'   \item Selects single predictors from \code{newdata} based on variables 
#'   listed in \code{object$keep}
#'   \item Adds combined predictors from paired variables in
#'    \code{object$pairs}
#' }
#' @export
.predict_multiDEGGs <- function(object,
                                newdata,
                                interaction.type = "ratio",
                                sep = ":",
                                ...) {
  result_data <- NULL
  keep <- object$keep
  
  # Add single predictors if present
  if(length(keep) > 0){
    single_preds <- newdata[, keep, drop = FALSE]
    result_data <- single_preds
  }
  
  # Add combined predictors if present
  if (nrow(object$pairs) > 0) {
    pairs <- object$pairs
    a <- newdata[, pairs[, 1], drop = FALSE]
    b <- newdata[, pairs[, 2], drop = FALSE]
    
    if (interaction.type == "ratio") {
      combined_preds <- a/b
    } else {
      combined_preds <- a*b
    }
    combined_preds[!is.finite(combined_preds)] <- 0
    colnames(combined_preds) <- paste(colnames(a), colnames(b), sep = sep)
    
    if (is.null(result_data)) {
      result_data <- combined_preds
    } else {
      result_data <- cbind(result_data, combined_preds)
    }
  }

  if (is.null(result_data)) {
    stop("The custom filter is empty")
  }
  return(result_data)
}


#' Wrapper of .predict_multiDEGGs for multiDEGGs_filter()
#'
#' This function generates predictions by creating a dataset with single 
#' and combined predictors based on the filtering results of a
#' multiDEGGs_filter model.
#'
#' @inherit .predict_multiDEGGs
#' @export
predict.multiDEGGs_filter <- function(object,
                                      newdata,
                                      interaction.type = "ratio",
                                      sep = ":", 
                                      ...) {
  .predict_multiDEGGs(object, newdata, interaction.type = "ratio",
                      sep = ":", ...)
}


#' Wrapper of .predict_multiDEGGs for multiDEGGs_filter_combined()
#'
#' This function generates predictions by creating a dataset with single 
#' and combined predictors based on the filtering results of a
#' multiDEGGs_filter model.
#'
#' @inherit .predict_multiDEGGs
#' @export
predict.multiDEGGs_filter_combined <- function(object,
                                               newdata,
                                               interaction.type = "ratio",
                                               sep = ":", 
                                               ...) {
  .predict_multiDEGGs(object, newdata, interaction.type = "ratio",
                      sep = ":", ...)
}


#' cat_parallel (from nestedcv)
#' 
#' Prints using shell echo from inside mclapply when run in Rstudio
#' @param ... to be passed to system()
cat_parallel <- function(...) {
  if (Sys.getenv("RSTUDIO") != "1") return()
  system(sprintf('echo "%s', paste0(..., '\\c"', collapse = "")))
}
