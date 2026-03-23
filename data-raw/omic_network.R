# The metapathway_gene_symbols object was obtained via the exportgraph() function of 
# the MITHrIL software (https://github.com/alaimos/mithril-standalone .
# This miRNA enriched network is then merged with all available networks in 
# Omnipath as follows: 
library(OmnipathR)

metapathway_gene_symbols <- readRDS("metapathway.RDS") 

all_dbs <- c(interaction_resources())
new_network <- omnipath_interactions(
  resources = all_dbs,
  organism = 9606
)
new_network <- new_network[!grepl("^COMPLEX", network$source) & 
                     !grepl("^COMPLEX", network$target), ]
new_network <- new_network[, c("source_genesymbol", "target_genesymbol")]
colnames(new_network) <- colnames(metapathway_gene_symbols)
new_network$link <- paste0(new_network$from, "_", new_network$to)
new_network <- as.data.frame(new_network)

metapathway_gene_symbols$link <- paste0(metapathway_gene_symbols$from, "_", metapathway_gene_symbols$to)

# merge
omic_network <- rbind(metapathway_gene_symbols, new_network)
omic_network <- subset(omic_network, !duplicated(omic_network$link))
omic_network$link <- NULL

tools::resaveRdaFiles("multiDEGGs/R/sysdata.rda", compress = "xz")

