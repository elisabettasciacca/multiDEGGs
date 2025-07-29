# multiDEGGs
## Differentially Expressed Gene-Gene pairs in multi omic data
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)
[![CRAN status](https://www.r-pkg.org/badges/version/multiDEGGs)](https://CRAN.R-project.org/package=multiDEGGs)
[![R-CMD-check](https://github.com/elisabettasciacca/multiDEGGs/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/elisabettasciacca/multiDEGGs/actions/workflows/R-CMD-check.yaml)  
[![downloads](https://cranlogs.r-pkg.org/badges/grand-total/multiDEGGs)](https://cran.r-project.org/package=multiDEGGs)
[![last commit](https://img.shields.io/github/last-commit/elisabettasciacca/multiDEGGs.svg)](https://github.com/elisabettasciacca/multiDEGGs/commits/master)

# multiDEGGs <a href="https://elisabettasciacca.github.io/multiDEGGs/"><img src="man/figures/logo.png" align="right" height="137" alt="multiDEGGs website" /></a>

The multiDEGGs package test for differential gene-gene correlations across different groups of samples in multi omic data.  
Specific gene-gene interactions can be explored and gene-gene pair regression plots can be interactively shown.   

### Installation 
Install from CRAN:    
`install.packages("multiDEGGs")`
  
Install from Github:  
`devtools::install_github("elisabettasciacca/multiDEGGs")`

### Example  
Load package and sample data   
```r
library(multiDEGGs)  
data("synthetic_metadata")  
data("synthetic_rnaseqData")  
data("synthetic_proteomicData")
data("synthetic_OlinkData")   
```
  
Generate differential networks:   
```r
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
                                 cores = 2)
```
  
Visualise interactively (will open a shiny interface)   
```r
View_diffNetworks(deggs_object)
```  
  
Get a table listing all the significant interactions found in each category  
```r
get_multiOmics_diffNetworks(deggs_object, sig_threshold = 0.05)
```
   
Plot differential regression fits for a single interaction  
`plot_regressions(deggs_object,
                 assayDataName = "RNAseq",
                 gene_A = "MTOR", 
                 gene_B = "AKT2",
                 legend_position = "bottomright")`
                  
## Citation
```r
citation("multiDEGGs")
```