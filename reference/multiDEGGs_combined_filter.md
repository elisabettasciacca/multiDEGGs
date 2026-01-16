# Combined multiDEGGs filter

This function can be passed to the `modifyX` parameter of
[nestcv.train](https://rdrr.io/pkg/nestedcv/man/nestcv.train.html) or
[nestcv.glmnet](https://rdrr.io/pkg/nestedcv/man/nestcv.glmnet.html) to
use one of the available statistical filters (t-test, wilcoxon, etc.) in
combination with multiDEGGs. Single predictors will be selected by the
selected statistical filter an paired predictors will be added by
multiDEGGs.

## Usage

``` r
multiDEGGs_combined_filter(
  y,
  x,
  filter_method = "ttest",
  nfilter,
  dynamic_nfilter = TRUE,
  keep_single_genes = FALSE,
  ...
)
```

## Arguments

- y:

  Numeric vector or factor. Response variable (outcome), i.e. the
  'metadata' named vector, as passed by
  [nestcv.train](https://rdrr.io/pkg/nestedcv/man/nestcv.train.html) or
  [nestcv.glmnet](https://rdrr.io/pkg/nestedcv/man/nestcv.glmnet.html).

- x:

  Predictor variables, i.e. the assayData matrix with genes in columns
  and IDs in rows, as passed by
  [nestcv.train](https://rdrr.io/pkg/nestedcv/man/nestcv.train.html) or
  [nestcv.glmnet](https://rdrr.io/pkg/nestedcv/man/nestcv.glmnet.html).

- filter_method:

  Character string. Statistical filtering method to be used in
  combination with multiDEGGs for sigle feature selection. Options are:
  "ttest", "wilcoxon", "ranger", "glmnet", "pls".

- nfilter:

  Integer. Maximum number of features to select.

- dynamic_nfilter:

  Logical. If `TRUE` `nfilter` will limit the number of features
  selected by the statistical filter and the feature space will be
  augmented by adding ALL the paired predictors found by multiDEGGs. If
  `FALSE` `nfilter` will limit the total number of predictors, with
  approximately half allocated to pairs and half to single genes.

- keep_single_genes:

  Logical. When `dynamic_nfilter = TRUE`, determines whether to include
  single genes selected by multiDEGGs (i.e. the single variables
  included in the differential pairs) in addition to those from the
  statistical filter. Default is FALSE.

- ...:

  Additional arguments passed to the filtering functions.

## Value

An object of class "multiDEGGs_filter" containing:

- keep:

  Character vector of selected single gene names

- pairs:

  Data frame of selected gene pairs with interaction information

## Details

The function operates in two modes:

**Dynamic Filtering (dynamic_nfilter = TRUE):**

- Selects `nfilter` single genes using the specified statistical method

- Finds all significant gene pairs using multiDEGGs

- Total predictors = nfilter single genes + number of significant pairs

  - If `keep_single_genes = TRUE`, also includes single genes obtained
    from pairs found by multiDEGGs

**Balanced Selection (dynamic_nfilter = FALSE):**

- Allocates approximately half of `nfilter` to gene pairs

- Remaining slots filled with single genes from the statistical filter

- If fewer pairs are found than allocated, compensates by selecting more
  single genes

- Ensures consistent total number of predictors across outer folds

The statistical filtering methods include:

- `"ttest"`: Two-sample t-test for differential expression

- `"wilcoxon"`: Wilcoxon rank-sum test

- `"ranger"`: Random Forest variable importance

- `"glmnet"`: Elastic net regularization

- `"pls"`: Partial Least Squares variable importance

## Examples

``` r
library(nestedcv)
data("synthetic_metadata")
data("synthetic_rnaseqData")

# fit a regularized linear model
# note that nfilter, n_outer_folds, n_inner_folds are set low to keep the
# example lightweight. Adjust these values as needed for your use case.
if (FALSE) { # \dontrun{
fit.glmnet <- nestedcv::nestcv.glmnet(
y = as.numeric(synthetic_metadata$response),
x =  t(synthetic_rnaseqData),
modifyX = "multiDEGGs_combined_filter",
modifyX_options = list(filter_method = "ttest", 
                       nfilter = 5,
                       dynamic_nfilter = TRUE, 
                       keep_single_genes = FALSE),
modifyX_useY = TRUE,
n_outer_folds = 4,
n_inner_folds = 4)

summary(fit.glmnet)
} # }

# fit a random forest model
# NOTE: nfilter, n_outer_folds, n_inner_folds are set low to keep the
# example lightweight. Adjust these values as needed for your use case.
fit.rf <- nestedcv::nestcv.train(
  y = synthetic_metadata$response,
  x = t(synthetic_rnaseqData),
  method = "rf",
  modifyX = "multiDEGGs_combined_filter",
  modifyX_options = list(filter_method = "ttest", 
                         nfilter = 5,
                         dynamic_nfilter = TRUE, 
                         keep_single_genes = FALSE),
  modifyX_useY = TRUE,
  n_outer_folds = 2,
  n_inner_folds = 2
)
#> Fitting final model using CV on whole data
#> Loading required package: ggplot2
#> Loading required package: lattice
#> Duration: 0.8327174 secs

fit.rf$summary
#>                Reference
#> Predicted       Non_responder Responder
#>   Non_responder            57         2
#>   Responder                 1        40
#> 
#>               AUC            Accuracy   Balanced accuracy   
#>            0.9963              0.9700              0.9676   
```
