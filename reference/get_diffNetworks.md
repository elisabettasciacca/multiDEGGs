# Generate multi-omic differential networks

Generate a multi-layer differential network with interaction p values

## Usage

``` r
get_diffNetworks(
  assayData,
  metadata,
  category_variable = NULL,
  regression_method = "lm",
  category_subset = NULL,
  network = NULL,
  percentile_vector = seq(0.35, 0.98, by = 0.05),
  padj_method = "bonferroni",
  show_progressBar = TRUE,
  verbose = TRUE,
  cores = parallel::detectCores()/2
)
```

## Arguments

- assayData:

  a matrix or data.frame (or list of matrices or data.frames for
  multi-omic analysis) containing normalised assay data. Sample IDs must
  be in columns and probe IDs (genes, proteins...) in rows. For multi
  omic analysis, it is highly recommended to use a named list of data.
  If unnamed, sequential names (assayData1, assayData2, etc.) will be
  assigned to identify each matrix or data.frame.

- metadata:

  a named vector, matrix, or data.frame containing sample annotations or
  categories. If matrix or data.frame, each row should correspond to a
  sample, with columns representing different sample characteristics
  (e.g., treatment group, condition, time point). The colname of the
  sample characteristic to be used for differential analysis must be
  specified in `category_variable`. Rownames must match the sample IDs
  used in assayData. If named vector, each element must correspond to a
  sample characteristic to be used for differential analysis, and names
  must match sample IDs used in the colnames of `assayData`. Continuous
  variables are not allowed.

- category_variable:

  when metadata is a matrix or data.frame this is the column name of
  `metadata` that contains the sample annotations to be used for
  differential analysis

- regression_method:

  whether to use robust linear modelling to calculate link p values.
  Options are 'lm' (default) or 'rlm'. The lm implementation is faster
  and lighter.

- category_subset:

  optional character vector indicating a subset of categories from the
  category variable. If not specified, all categories in
  `category_variable` will be used.

- network:

  network of biological interactions provided by the user. The network
  must be provided in the form of a table of class data.frame with only
  two columns named "from" and "to". If NULL (default) a network of
  10,537 molecular interactions obtained from KEGG, mirTARbase,
  miRecords and transmiR will be used. This has been obtained via the
  `exportgraph` function of the MITHrIL tool (Alaimo et al., 2016).

- percentile_vector:

  a numeric vector specifying the percentiles to be used in the
  percolation analysis. By default, it is defined as
  `seq(0.35, 0.98, by = 0.05)`, which generates a sequence of
  percentiles starting at 0.35, meaning that targets (genes/proteins...)
  whose expression value is under the 35th percentile of the whole
  matrix will be excluded. This threshold can be modified by specifying
  a different starting point for `seq`. For a more granular percolation
  analysis an higher optimisation of the algorithm, `by = 0.05` can be
  modified in favour of lower values, but this will increase the
  computational time.

- padj_method:

  a character string indicating the p values correction method for
  multiple test adjustment. It can be either one of the methods provided
  by the `p.adjust` function from `stats` (bonferroni, BH, hochberg,
  etc.) or "q.value" for Storey's q values, or "none" for unadjusted p
  values. When using "q.value" the `qvalue` package must be installed
  first.

- show_progressBar:

  logical. Whether to display a progress bar during execution. Default
  is TRUE.

- verbose:

  logical. Whether to print detailed output messages during processing.
  Default is TRUE

- cores:

  number of cores to use for parallelisation.

## Value

a `deggs` object containing differential networks incorporating p values
or adjusted p values for each link.

## Examples

``` r
# Single omic analysis:
data("synthetic_metadata")
data("synthetic_rnaseqData")
deggs_object <- get_diffNetworks(assayData = synthetic_rnaseqData,
                                 metadata = synthetic_metadata,
                                 category_variable = "response",
                                 regression_method = "lm",
                                 padj_method = "bonferroni",
                                 verbose = FALSE,
                                 show_progressBar = FALSE,
                                 cores = 1)
 
# multi-omic analysis: 
data("synthetic_metadata")
data("synthetic_rnaseqData")
data("synthetic_proteomicData")
data("synthetic_OlinkData")
assayData_list <- list("RNAseq" = synthetic_rnaseqData,
                       "Proteomics" = synthetic_proteomicData,
                       "Olink" = synthetic_OlinkData)
deggs_object <- get_diffNetworks(assayData = assayData_list,
                                 metadata = synthetic_metadata,
                                 category_variable = "response",
                                 regression_method = "lm",
                                 padj_method = "bonferroni",
                                 verbose = FALSE,
                                 show_progressBar = FALSE,
                                 cores = 1)
                                 
# to use only certain categories for comparison: 
# let's randomly add another level of response to the example metadata
synthetic_metadata$response <- as.character(synthetic_metadata$response)
indices <- sample(1:nrow(synthetic_metadata), 20, replace = FALSE) 
synthetic_metadata$response[indices] <- "Moderate response"
deggs_object <- get_diffNetworks(assayData = assayData_list,
                                 metadata = synthetic_metadata,
                                 category_variable = "response",
                                 category_subset = c("Responder", 
                                                     "Non_responder"),
                                 regression_method = "lm",
                                 verbose = FALSE,
                                 show_progressBar = FALSE,
                                 cores = 1)
                                 
# to be more generous on the targets to be excluded, and lower the expression 
# level threshold to the 25th percentile (or lower): 
deggs_object <- get_diffNetworks(assayData = assayData_list,
                                 metadata = synthetic_metadata,
                                 category_variable = "response",
                                 category_subset = c("Responder", 
                                                     "Non_responder"),
                                 regression_method = "lm",
                                 percentile_vector = seq(0.25, 0.98, by = 0.05),
                                 verbose = FALSE,
                                 show_progressBar = FALSE,
                                 cores = 1)                                  
```
