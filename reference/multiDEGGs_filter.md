# multiDEGGs_filter

Function to be passed to the `modifyX` parameter of
[nestcv.train](https://rdrr.io/pkg/nestedcv/man/nestcv.train.html) or
[nestcv.glmnet](https://rdrr.io/pkg/nestedcv/man/nestcv.glmnet.html) to
allow nested feature selection and augmentation via differential network
analysis with multiDEGGs.

## Usage

``` r
multiDEGGs_filter(y, x, keep_single_genes = FALSE, nfilter = NULL)
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

- keep_single_genes:

  Logical, default FALSE. If TRUE, the function will return unique
  individual genes along with significant pairs.

- nfilter:

  Integer. Maximum total number of predictors to return. When
  `keep_single_genes = TRUE`, this parameter limits the combined count
  of unique and paired predictors (i.e., length(keep_DEGGs) +
  nrow(pairs) \<= nfilter). Predictors are included from most to least
  significant: for each pair, both the pair itself and the new unique
  variables are included until the `nfilter` threshold is reached. When
  `keep_single_genes = FALSE`, `nfilter` only limits the number of pairs
  returned. If `NULL`, no filtering is applied and the total number of
  predictors will depend on how many significantly different
  interactions are detected by multiDEGGs in that fold.

## Value

a list containing two types of predictors:

- single predictors (stored in the 'keep_DEGGs' vector)

- paired predictors (stored in the 'pairs' dataframe) Note that
  `nfilter` limits the maximum number of engineered features returned,
  however this number might be lower and will depend on how many
  significantly different interactions will be found in each fold by
  multiDEGGs. **If no significantly different interactions are found the
  function will print a '0' and switch to t-test for that fold.**

## Examples

``` r
library(nestedcv)
data("synthetic_metadata")
data("synthetic_rnaseqData")

# fit a regularized linear model
# Note that nfilter, n_outer_folds, n_inner_folds are set low to keep the
# example lightweight. Adjust these values as needed for your use case.
if (FALSE) { # \dontrun{
fit.glmnet <- nestcv.glmnet(
  y = as.numeric(synthetic_metadata$response),
  x =  t(synthetic_rnaseqData),
  modifyX = "multiDEGGs_filter",
  modifyX_options = list(keep_single_genes = FALSE,
                         nfilter = 5),
  modifyX_useY = TRUE,
  n_outer_folds = 4,
  n_inner_folds = 4)

summary(fit.glmnet)
} # }

# fit a random forest model:
# note that nfilter, n_outer_folds, n_inner_folds are set low to keep the
# example lightweight. Adjust these values as needed for your use case.
fit.rf <- nestcv.train(
  y = synthetic_metadata$response,
  x = t(synthetic_rnaseqData),
  method = "rf",
  modifyX = "multiDEGGs_filter",
  modifyX_options = list(keep_single_genes = FALSE,
                         nfilter = 5),
  modifyX_useY = TRUE,
  n_outer_folds = 2,
  n_inner_folds = 2
)
#> Fitting final model using CV on whole data
#> Duration: 0.7319274 secs

fit.rf$summary
#>                Reference
#> Predicted       Non_responder Responder
#>   Non_responder            54        12
#>   Responder                 4        30
#> 
#>               AUC            Accuracy   Balanced accuracy   
#>            0.8996              0.8400              0.8227   
```
