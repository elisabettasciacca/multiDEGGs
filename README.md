# multiDEGGs
## Differentially Expressed Gene-Gene pairs in multi omic data
The multiDEGGs package test for differential gene-gene correlations across different groups of samples in multi omic data.  
Specific gene-gene interactions can be explored and gene-gene pair regression plots can be interactively shown.   

### Installation 
Install from CRAN:    
`install.packages("multiDEGGs")`
  
Install from Github:  
`devtools::install_github("elisabettasciacca/multiDEGGs")`

### Example  
Load package and sample data   
`library(multiDEGGs)
data("synthetic_metadata")
data("synthetic_rnaseqData")
data("synthetic_proteomicData")
data("synthetic_OlinkData")`  

Generate differential networks   
`assayData_list <- list("RNAseq" = synthetic_rnaseqData,
                       "Proteomics" = synthetic_proteomicData,
                       "Olink" = synthetic_OlinkData)

deggs_object <- get_diffNetworks(assayData = assayData_list,
                                 metadata = synthetic_metadata,
                                 category_variable = "response",
                                 regression_method = "lm",
                                 padj_method = "bonferroni",
                                 verbose = FALSE,
                                 show_progressBar = FALSE,
                                 cores = 2)` 
  
Visualise interactively (will open a shiny interface)   
`View_diffNetworks(deggs_object)`  
  
Get a table listing all the significant interactions found in each category  
`get_multiOmics_diffNetworks(deggs_object, sig_threshold = 0.05)`  
   
Plot differential regression fits for a single interaction  
`plot_regressions(deggs_object,
                 assayDataName = "RNAseq",
                 gene_A = "MTOR", 
                 gene_B = "AKT2",
                 legend_position = "bottomright")`
                  
## Citation
```{r}
citation("multiDEGGs")
```