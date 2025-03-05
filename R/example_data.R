#' TCGA Breast Invasive Carcinoma sample data
#'
#' A dataset containing sample data for 105 breast tissue biopsies from the
#' TCGA-BRCA project
#'
#' @format A data frame with 105 rows and 1 variable (IDs are in rownames):
#' \describe{
#'   \item{SUBTYPE}{The breast cancer subtype}
#' }
"BRCA_metadata"

#' TCGA Breast Invasive Carcinoma (BRCA) normalised gene expression data
#'
#' A dataset containing the z scored gene expression data for 105 breast
#' biopsies from the TCGA-BRCA project
#'
#' @format A data frame with 19737 rows representing normalised gene counts
#' and 105 columns representing samples.
"BRCA_normCounts"
