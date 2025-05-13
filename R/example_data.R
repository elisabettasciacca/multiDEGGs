#' Synthetic clinical data
#'
#' A dataset containing sample clinical data for 100 patients with 40% 
#' response rate
#' 
#' @format A data frame with 100 rows and 4 columns (IDs are in rownames):
#' \describe{
#'   \item{response}{The response outcome, to be used for differential analysis}
#' }
"synthetic_metadata"

#' Synthetic RNA-seq count data
#'
#' Synthetic RNA-seq data after log2 normalisation
#'
#' @format A data frame with xx rows (genes) xx columns (patients IDs, matching 
#' the metadata rownames).
"synthetic_rnaseqData"

#' Synthetic RNA-seq count data
#'
#' Synthetic RNA-seq data after log2 normalisation
#'
#' @format A data frame with xx rows (proteins) xx columns (patients IDs).
"synthetic_proteomicData"

#' Synthetic RNA-seq count data
#'
#' Synthetic RNA-seq data after log2 normalisation
#'
#' @format A data frame with xx rows (proteins) xx columns (patients IDs).
"synthetic_OlinkData"
