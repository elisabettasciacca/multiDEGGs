rm(list = ls())
source("Percolation.R")
library("dplyr")
library("stringr")

# removing self loops from the two internal networks ("metapathways")
metapathway_gene_symbols <- metapathway_gene_symbols[which(!(
  metapathway_gene_symbols$from == metapathway_gene_symbols$to)), ]

metapathway_entrez_IDs <- metapathway_entrez_IDs[which(!(
  metapathway_entrez_IDs$from == metapathway_entrez_IDs$to)), ]
#save to rda file
save(metapathway_gene_symbols, metapathway_entrez_IDs, file = "sysdata.rda")

# import and refactor of normalised_counts dataset
normalised_counts <- read.table(
  file = 'data_RNA_Seq_v2_mRNA_median_Zscores_normals.txt',
  na.strings = c("", "NA"),
  sep = '\t', header = TRUE
)[, -1]
rownames(normalised_counts) <- normalised_counts[, 1]
normalised_counts[, 1] <- NULL

# renaming normalised_counts column names
normalised_counts <- normalised_counts %>%
  rename_all(funs(str_replace_all(., '\\.11', '') %>%
    str_replace_all(., '\\.', '-')))

# removing null genes
normalised_counts <- na.omit(normalised_counts)

# import and refactor of metadata dataset
metadata <- read.table(
  file = 'data_clinical_patient.tsv',
  na.strings = c("", "NA"), sep = '\t',
  header = TRUE
)[, -c(3:38)]
rownames(metadata) <- metadata[, 1]
metadata[, 1] <- NULL

# removing null samples
print(nrow(na.omit(metadata)))
metadata <- na.omit(metadata)

# align data
metadata <- metadata[colnames(normalised_counts), , drop = FALSE]
metadata <- na.omit(metadata)
normalised_counts <- normalised_counts[, rownames(metadata)]


usethis::use_data(THCA_metadata, THCA_norm_counts)
