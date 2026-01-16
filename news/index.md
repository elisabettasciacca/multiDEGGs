# Changelog

## multiDEGGs 1.1.2

CRAN release: 2026-01-15

## multiDEGGs 1.1.2

CRAN release: 2026-01-15

Patch:
[`get_diffNetworks_singleOmic()`](https://elisabettasciacca.github.io/multiDEGGs/reference/get_diffNetworks_singleOmic.md)
now double checks that `assayData` and `metadata` are aligned for the
layer

## multiDEGGs 1.1.1

CRAN release: 2025-10-24

Minor fixes in documentation

## multiDEGGs 1.1.0

CRAN release: 2025-07-29

#### New features for feature augmentation in ML

Two new functions are provided for nested feature engineering. To use
them in combination with the `nestedcv` package their name must be
passed to the `modifyX` parameter of
[`nestcv.glmnet()`](https://rdrr.io/pkg/nestedcv/man/nestcv.glmnet.html)
or
[`nestcv.train()`](https://rdrr.io/pkg/nestedcv/man/nestcv.train.html).

- The
  [`multiDEGGs_filter()`](https://elisabettasciacca.github.io/multiDEGGs/reference/multiDEGGs_filter.md)
  function performs feature selection based entirely on differential
  network analysis.
- The
  [`multiDEGGs_combined_filter()`](https://elisabettasciacca.github.io/multiDEGGs/reference/multiDEGGs_combined_filter.md)
  function combines traditional statistical feature selection (5
  options) with differential network analysis.
- Internally the two
  [`predict.multiDEGGs_filter()`](https://elisabettasciacca.github.io/multiDEGGs/reference/predict.multiDEGGs_filter.md)
  and `predict.multiDEGGs_combined_filter()` S3 methods generate
  predictions by creating a dataset with single and combined predictors
  based on the filtering results of a `multiDEGGs_filter` model.
- The vignette has been updated to showcase the new feature

## multiDEGGs 1.0.0

CRAN release: 2025-06-05

#### Initial Release

- First public version of `multiDEGGs`
- Provides tools for differential network analysis.
- Can be easily integrated in machine learning pipelines as feature
  selection method.
- Supports both single and multi omic analyses.
- Compatible with R \>= 4.4.
