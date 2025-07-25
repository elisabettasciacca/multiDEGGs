# multiDEGGs_filter - test -------------------------------------------------
test_that("test multiDEGGs_filter", {
  skip_if_not_installed("nestedcv")
  set.seed(123)
  data("synthetic_metadata")
  data("synthetic_rnaseqData")

  res <- multiDEGGs_filter(y = as.numeric(synthetic_metadata$response),
                    x =  t(synthetic_rnaseqData))
  
  expect_s3_class(res, "multiDEGGs_filter")
  expect_equal(ncol(res$pairs), 4)
  expect_equal(nrow(res$pairs), 7)
  expect_true("IL1B" %in% res$pairs$from)
  expect_true("IL1R2" %in% res$pairs$to)
  expect_true("FASLG" %in% res$pairs$from)
  expect_true("FAS" %in% res$pairs$to)
  
  fit.glmnet <- nestedcv::nestcv.glmnet(
    y = as.numeric(synthetic_metadata$response),
    x =  t(synthetic_rnaseqData),
    modifyX = "multiDEGGs_filter",
    modifyX_options = list(keep_single_genes = FALSE,
                           nfilter = 20),
    modifyX_useY = TRUE,
    n_outer_folds = 5,
    n_inner_folds = 6)
  
  fit.glmnet.summary <- summary(fit.glmnet)
  
  expect_true("TNF:TNFRSF1A" %in% names(fit.glmnet.summary[["coef"]]))
  
  fit.svm <- nestedcv::nestcv.train(
    y = synthetic_metadata$response,
    x = t(synthetic_rnaseqData),
    method = "svmRadial",
    modifyX = "multiDEGGs_filter",
    modifyX_options = list(keep_single_genes = TRUE,
                           nfilter = 6),
    modifyX_useY = TRUE,
    n_outer_folds = 5,
    n_inner_folds = 6,
    verbose = FALSE
  )
  
  expect_s3_class(fit.svm, "nestcv.train")
  expect_equal(length(fit.svm[["final_vars"]]), 6)
})

# multiDEGGs_combined_filter - test -------------------------------------------------
test_that("test multiDEGGs_combined_filter", {
  skip_if_not_installed("nestedcv")
  set.seed(123)
  data("synthetic_metadata")
  data("synthetic_rnaseqData")
  
  res <- multiDEGGs_combined_filter(y = as.numeric(synthetic_metadata$response),
                           x =  t(synthetic_rnaseqData), 
                           filter_method = "ttest",
                           nfilter = 10,
                           dynamic_nfilter = FALSE)
  
  expect_s3_class(res, "multiDEGGs_filter")
  expect_true("MAPK3" %in% res[["keep"]])
  expect_true("TGFB3" %in% res[["keep"]])
  expect_equal(ncol(res$pairs), 4)
  expect_equal(nrow(res$pairs), 5)
  expect_true("FANCD2" %in% res$pairs$from)
  expect_true("FAN1" %in% res$pairs$to)
  
  res2 <- multiDEGGs_combined_filter(y = as.numeric(synthetic_metadata$response),
                                    x =  t(synthetic_rnaseqData), 
                                    filter_method = "wilcoxon",
                                    nfilter = 10,
                                    dynamic_nfilter = TRUE)
  
  
  expect_s3_class(res2, "multiDEGGs_filter")
  expect_equal(length(res2[["keep"]]), 10)
  expect_true(all(c("TNF", "MAP2K2", "MAPK3", "TNFRSF1A", "MTOR", "FAN1",
                    "IL1B") %in% res2[["keep"]]))
  expect_equal(ncol(res2$pairs), 4)
  expect_equal(nrow(res2$pairs), 7)
  expect_true("IL1B" %in% res2$pairs$from)
  expect_true("IL1R2" %in% res2$pairs$to)
  
  res3 <- multiDEGGs_combined_filter(y = as.numeric(synthetic_metadata$response),
                                     x =  t(synthetic_rnaseqData), 
                                     filter_method = "pls",
                                     nfilter = 10,
                                     dynamic_nfilter = TRUE,
                                     keep_single_genes = TRUE)
  
  expect_s3_class(res3, "multiDEGGs_filter")
  expect_equal(length(res3[["keep"]]), 10)
  expect_true(all(c("MAP2K2", "TNFRSF1A", "MTOR", "FANCD2", "AKT2", 
                    "FAN1" ) %in% res3[["keep"]]))
  expect_equal(ncol(res3$pairs), 4)
  expect_equal(nrow(res3$pairs), 7)
  expect_true("IL1B" %in% res3$pairs$from)
  expect_true("TNFRSF1A" %in% res3$pairs$to)
  
  fit.glmnet <- nestedcv::nestcv.glmnet(
    y = as.numeric(synthetic_metadata$response),
    x =  t(synthetic_rnaseqData),
    modifyX = "multiDEGGs_combined_filter",
    modifyX_options = list(filter_method = "ttest", 
                           nfilter = 20,
                           dynamic_nfilter = TRUE, 
                           keep_single_genes = FALSE),
    modifyX_useY = TRUE,
    n_outer_folds = 5,
    n_inner_folds = 6)
  
  fit.glmnet.summary <- summary(fit.glmnet)
  
  expect_equal(length(fit.glmnet[["final_vars"]]), 21)
  expect_equal(nrow(fit.glmnet[["final_coef"]]), 15)
  
  expect_true("MAP2K2:MAPK3" %in% fit.glmnet[["final_vars"]])
  
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
    n_outer_folds = 5,
    n_inner_folds = 6,
    verbose = FALSE
  )
  
  expect_s3_class(fit.rf, "nestcv.train")
  expect_equal(fit.rf[["summary"]][["metrics"]][["AUC"]], 0.99, 
               tolerance = 0.01)
  expect_equal(fit.rf[["summary"]][["metrics"]][["Accuracy"]], 0.98, 
               tolerance = 0.01)
  expect_equal(length(fit.rf[["final_vars"]]), 12)
})