# get_diffNetworks - test -------------------------------------------------
test_that("test get_diffNetworks multi_omics", {
  # let's prepare a multi-omic example
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
  
  # check the structure of the result is correct 
  expect_s3_class(deggs_object, "deggs")
  expect_named(deggs_object, c("diffNetworks", "assayData", "metadata", 
                               "regression_method", "category_subset",  
                               "padj_method"))
  expect_equal(deggs_object$assayData, assayData_list)
  expect_equal(length(deggs_object$assayData), length(assayData_list))
  
  expect_s3_class(deggs_object$metadata, "factor")
  expect_type(deggs_object$regression_method, "character")
  expect_type(deggs_object$padj_method, "character")
  
  # check the content hasn't changed 
  expect_equal(deggs_object[["diffNetworks"]][["RNAseq"]][["sig_pvalues_count"]],
               7)
  expect_equal(nrow(deggs_object[["diffNetworks"]][["RNAseq"]][["Responder"]]),
               3)
  expect_equal(deggs_object[["diffNetworks"]][["Proteomics"]][["Responder"]],
               "No specific links for this category.")
  expect_equal(ncol(deggs_object[["diffNetworks"]][["Olink"]][["Non_responder"]]),
               4)
  expect_equal(deggs_object[["diffNetworks"]][["Olink"]][["Responder"]],
               "No specific links for this category.")
  
  # even more granular check on the differential links found 
  expect_true("TNF" %in% 
        deggs_object[["diffNetworks"]][["RNAseq"]][["Non_responder"]][["from"]])
  expect_true("TNFRSF1A" %in% 
         deggs_object[["diffNetworks"]][["RNAseq"]][["Non_responder"]][["to"]])
})


# get_multiOmics_diffNetworks - test -------------------------------------------------
test_that("test get_multiOmics_diffNetworks multi_omics", {
# let's prepare a multi- omic example
  data("synthetic_metadata")
  data("synthetic_rnaseqData")
  data("synthetic_OlinkData")
  
  assayData_list <- list("RNAseq" = synthetic_rnaseqData,
                         "Olink" = synthetic_OlinkData)
  
  deggs_object <- get_diffNetworks(assayData = assayData_list,
                                   metadata = synthetic_metadata,
                                   category_variable = "response",
                                   regression_method = "rlm",
                                   padj_method = "q.value",
                                   verbose = FALSE,
                                   show_progressBar = FALSE,
                                   cores = 1)
res <- get_multiOmics_diffNetworks(deggs_object)

# check the structure of the result is correct 
expect_type(res, "list")

# check the content hasn't changed 
expect_equal(length(res), 2)
expect_equal(ncol(res$`Non_responder`), 5)
expect_equal(colnames(res$Responder), c("from", "to", "p.value", "p.adj",
                                        "layer"))
expect_equal(nrow(res$Responder), 3)

# even more granular check on the differential links found 
expect_true("AKT2" %in% res$`Non_responder`$from)
expect_true("MTOR" %in% res$`Non_responder`$to)

expect_true("FANCD2" %in% res$Responder$from)
expect_true("FAN1" %in% res$Responder$to)

})


# get_diffNetworks - test -------------------------------------------------
test_that("test get_diffNetworks single_omics", {
  # let's prepare a single omic example
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
  res <- get_sig_deggs(deggs_object)
  
  # check the structure of the result is correct 
  expect_s3_class(res, "data.frame")
  
  # check the content hasn't changed 
  expect_equal(ncol(res), 4)
  expect_equal(nrow(res), 7)
  
  # even more granular check on the differential links found 
  expect_true("TNF" %in% res$from)
  expect_true("TNFRSF1A" %in% res$to)
  
  expect_true("TGFB3" %in% res$from)
  expect_true("TGFBR1" %in% res$to)
  
  expect_true("IL1B" %in% res$from)
  expect_true("IL1R2" %in% res$to)
  
})
