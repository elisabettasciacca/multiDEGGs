test_that("test get_multiOmics_diffNetworks multi_omics", {
  # let's prepare a single omic example
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

  expect_silent(plot_regressions(deggs_object,
                   assayDataName = "RNAseq",
                   gene_A = "MTOR", 
                   gene_B = "AKT2",
                   legend_position = "bottomright"))
  
  expect_type(plot_regressions(deggs_object,
                                 assayDataName = "RNAseq",
                                 gene_A = "MTOR", 
                                 gene_B = "AKT2",
                                 legend_position = "bottomright"), "list")
  
  expect_no_warning(plot_regressions(deggs_object,
                                 assayDataName = "RNAseq",
                                 gene_A = "MTOR", 
                                 gene_B = "AKT2",
                                 legend_position = "bottomright"))
})
  
test_that("plot on graphical device", {
  skip_on_cran()  
  skip_on_ci()   
  skip_if_not(interactive())  
  
  expect_true({
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
    plot_regressions(deggs_object,
                     assayDataName = "RNAseq",
                     gene_A = "MTOR", 
                     gene_B = "AKT2",
                     legend_position = "bottomright")
    length(dev.list()) >= 1
  })
})

test_that("test View_diffNetworks multi_omics", {
  # let's prepare a multi-omic example
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
  
  res <- View_diffNetworks(deggs_object)
  
  expect_s3_class(res, "shiny.appobj")
  expect_type(res, "list")
  expect_true(all(c("httpHandler", "serverFuncSource", "onStart", "options") %in% 
                names(res)))
  
})