# multiDEGGs 1.1.0
### New features for feature augmentation in ML
Two new functions are provided for nested feature engineering. To use them in 
combination with the `nestedcv` package their name must be passed to the 
`modifyX` parameter of `nestcv.glmnet()` or `nestcv.train()`. 
  
- The `multiDEGGs_filter()` function performs feature selection based entirely 
on differential network analysis. 
- The `multiDEGGs_combined_filter()` function combines traditional statistical
feature selection (5 options) with differential network analysis. 
- Internally the two `predict.multiDEGGs_filter()` and 
`predict.multiDEGGs_combined_filter()` S3 methods generate predictions by 
creating a dataset with single and combined predictors based on the filtering 
results of a `multiDEGGs_filter` model.
- The vignette has been updated to showcase the new feature

# multiDEGGs 1.0.0
### Initial Release
- First public version of `multiDEGGs`
- Provides tools for differential network analysis.
- Can be easily integrated in machine learning pipelines as feature selection method.
- Supports both single and multi omic analyses.
- Compatible with R >= 4.4.
