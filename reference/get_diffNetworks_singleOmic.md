# Generate differential networks for single omic analysis

Generate differential networks for single omic analysis

## Usage

``` r
get_diffNetworks_singleOmic(
  assayData,
  assayDataName,
  metadata,
  regression_method,
  network,
  percentile_vector,
  padj_method,
  show_progressBar,
  verbose,
  cores
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

- assayDataName:

  name of the assayData, to identify which omic is.

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

- regression_method:

  whether to use robust linear modelling to calculate link p values.
  Options are 'lm' (default) or 'rlm'. The lm implementation is faster
  and lighter.

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

a list of differential networks, one per category
