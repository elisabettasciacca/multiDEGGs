# Calculate the p values for specific category network samples

Calculate the p values for specific category network samples

## Usage

``` r
calc_pvalues_network(
  assayData,
  metadata,
  padj_method,
  categories_length,
  regression_method = "lm",
  category_network
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

- padj_method:

  a character string indicating the p values correction method for
  multiple test adjustment. It can be either one of the methods provided
  by the `p.adjust` function from `stats` (bonferroni, BH, hochberg,
  etc.) or "q.value" for Storey's q values, or "none" for unadjusted p
  values. When using "q.value" the `qvalue` package must be installed
  first.

- categories_length:

  integer number indicating the number of categories

- regression_method:

  whether to use robust linear modelling to calculate link p values.
  Options are 'lm' (default) or 'rlm'. The lm implementation is faster
  and lighter.

- category_network:

  network table for a specific category

## Value

a list of p values
