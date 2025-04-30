# DEGGs
## Differentially Expressed Gene-Gene pairs
The multiDEGGs package test for differential gene-gene correlations across different groups of samples in multi omic data.  
Specific gene-gene interactions can be explored and gene-gene pair regression plots can be interactively shown.   

### Installation instructions 
To install from Github please use the following on your R console  
`devtools::install_github("elisabettasciacca/multiDEGGs", build_vignettes = TRUE)`

### Example  
Load package and sample data   
`library(multiDEGGs)  
data("BRCA_metadata")  
data("BRCA_normCounts")`  

Generate specific gene-gene networks for each subtype  
`deggs_object <- get_diffNetworks(assayData = syn.BRCA_normCounts,
                                       metadata = BRCA_metadata,
                                       category_variable = "SUBTYPE", 
                                       category_subset = c("BRCA_Her2", "BRCA_LumA"), 
                                       use_qvalues = TRUE, 
                                       verbose = TRUE,
                                       cores = 2)` 
  
Visualise  
`View_diffNetworks(deggs_object)`  
  
Get a table listing all the significant gene-gene interactions found in each subtype  
`get_sig_deggs(deggs_object)`  
   
Print differential regression fits for a single gene-gene interaction through the `print_regressions` function  
`print_regressions(gene_A = "NOTCH2", gene_B = "DTX4",
                  deggs_object = deggs_object,
                  legend_position = "bottomright")`
                  
## Citation

multiDEGGs was developed by Elisabetta Sciacca and supported by the bioinformatics team at 
[Experimental Medicine & Rheumatology department](https://www.qmul.ac.uk/whri/emr/) 
and [Centre for Translational Bioinformatics](https://www.qmul.ac.uk/c4tb/) (Queen Mary University London), in joint collaboration with the 
[Department of Clinical and Experimental Medicine at University of Catania](https://www.medclin.unict.it/en). 


If you use this package please cite as: 

```{r}
citation("multiDEGGs")
```

or:

> Sciacca, Elisabetta, et al. "DEGGs: an R package with shiny app for the identification of differentially expressed gene–gene interactions in high-throughput sequencing data." Bioinformatics 39.4 (2023): btad192.
