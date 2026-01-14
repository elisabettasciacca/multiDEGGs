# Predict method for multiDEGGs_filter objects

This function generates predictions by creating a dataset with single
and combined predictors based on the filtering results of a
multiDEGGs_filter model.

## Usage

``` r
.predict_multiDEGGs(
  object,
  newdata,
  interaction.type = "ratio",
  sep = ":",
  ...
)
```

## Arguments

- object:

  A fitted object of class `multiDEGGs_filter` containing filtering
  results with:

  keep

  :   Character vector of variable names to keep as single predictors

  pairs

  :   Data frame or matrix with two columns specifying pairs of
      variables to combine

- newdata:

  A data frame containing the new data for prediction. Must contain all
  variables specified in `object$keep` and `object$pairs`.

- interaction.type:

  Character string specifying how to combine the paired predictors.
  Options are:

  "ratio"

  :   Combine paired predictors by dividing the first variable by the
      second (a/b)

  other

  :   Combine paired predictors by multiplying the variables (a\*b)

  Default is "ratio".

- sep:

  Character string used as separator when creating column names for
  combined predictors. Default is ":".

- ...:

  Additional arguments passed to the generic function.

## Value

A data frame containing:

- Single predictors (if any are specified in `object$keep`)

- Combined predictors based on variable pairs and interaction type

## Details

The function processes the filtering results in two steps:

1.  Selects single predictors from `newdata` based on variables listed
    in `object$keep`

2.  Adds combined predictors from paired variables in `object$pairs`
