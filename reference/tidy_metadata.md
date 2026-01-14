# Tidying up of metadata. Samples belonging to undesidered categories (if specified) will be removed as well as categories with less than five samples, and NAs.

Tidying up of metadata. Samples belonging to undesidered categories (if
specified) will be removed as well as categories with less than five
samples, and NAs.

## Usage

``` r
tidy_metadata(
  category_subset = NULL,
  metadata,
  category_variable = NULL,
  verbose = FALSE
)
```

## Arguments

- category_subset:

  optional character vector indicating which categories are used for
  comparison. If not specified, all categories in `category_variable`
  will be used.

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

  column name in `metadata` (if data.frame or matrix) or NULL if
  `metadata` is already a named vector containing category information.

- verbose:

  Logical. Whether to print detailed output messages during processing.
  Default is FALSE.

## Value

a tidy named factor vector of sample annotations.
