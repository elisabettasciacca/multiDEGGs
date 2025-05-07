# Synthetic Data Generator for multi-omic Differential Network Analysis
# This script generates fake but realistic multi-omic datasets with:
# 1. Metadata with responders/non-responders
# 2. RNA-seq count matrix 
# 3. Proteomics matrix 
# 4. Olink matrix 
# All data is supposed to be log2 normalised

rm(list = ls())
set.seed(122)

# Define parameters ----------------------------------
n_patients <- 100  
n_genes <- 180      
n_proteins <- 130   
n_olink <- 92 
response_rate <- 0.4  
# Add some realistic clinical variables
age_range <- c(35, 75)
gender_ratio <- 0.55  # proportion female

# Use the metapathway gene network included in the packege ------------------
load("R/sysdata.rda")

# Extract unique genes from the metapathway
all_genes <- unique(c(metapathway_gene_symbols$from, metapathway_gene_symbols$to))

metadata_generation <- function() {
  patient_ids <- paste0("PT", sprintf("%03d", 1:n_patients))
  
  response <- sample(c("Responder", "Non-responder"), 
                     size = n_patients, 
                     replace = TRUE, 
                     prob = c(response_rate, 1-response_rate))
  
  # Add age with slight difference between groups
  # Responders tend to be slightly younger
  age <- numeric(n_patients)
  for (i in 1:n_patients) {
    if (response[i] == "Responder") {
      age[i] <- round(runif(1, age_range[1], age_range[2] - 5))
    } else {
      age[i] <- round(runif(1, age_range[1] + 3, age_range[2]))
    }
  }
  
  # Generate gender with slight imbalance
  gender <- sample(c("Female", "Male"), 
                   size = n_patients, 
                   replace = TRUE, 
                   prob = c(gender_ratio, 1-gender_ratio))
  
  # Create a data frame
  metadata <- data.frame(
    patient_id = patient_ids,
    response = response,
    age = age,
    gender = gender,
    stringsAsFactors = FALSE
  )
  
  return(metadata)
}

data_correlated_generation <- function(n, mean_x, mean_y, sd_x, sd_y, correlation, direction, rand_number) {
  # Apply the random number based on the direction passed from the parent function
  mean_x <- mean_x + (direction * rand_number)
  mean_y <- mean_y + (direction * rand_number)
  
  # Create covariance matrix
  sigma <- matrix(c(sd_x^2, sd_x*sd_y*correlation, sd_x*sd_y*correlation, sd_y^2), 2, 2)
  
  # Generate multivariate normal data
  data <- MASS::mvrnorm(n, mu = c(mean_x, mean_y), Sigma = sigma)
  return(data)
}

correlation_differential_generation <- function(responder_indices, nonresponder_indices, 
                                                mean_val, sd_val, 
                                                entity1_idx, entity2_idx,
                                                matrix_data) {
  
  # Decide on correlation pattern type
  pattern_type <- sample(1:4, 1)
  
  # Generate random number to be added/subtracted to mean values of responders 
  # and non-responders (between 1 and 3)
  rand_number <- runif(1, 1, 3)
  
  # Randomly decide direction (add to responders or non-responders)
  # 1 means add to responders, subtract from non-responders
  # -1 means subtract from responders, add to non-responders
  direction <- sample(c(1, -1), 1)
  
  if (pattern_type == 1) {
    # Pattern 1: Positive correlation in responders, negative in non-responders
    
    # For responders: strong positive correlation
    resp_data <- data_correlated_generation(
      n = length(responder_indices),
      mean_x = mean_val, mean_y = mean_val, 
      sd_x = sd_val, sd_y = sd_val,
      correlation = runif(1, 0.7, 0.9),  # Strong positive correlation
      direction = direction,
      rand_number = rand_number
    )
    
    # For non-responders: negative correlation
    nonresp_data <- data_correlated_generation(
      n = length(nonresponder_indices),
      mean_x = mean_val, mean_y = mean_val, 
      sd_x = sd_val, sd_y = sd_val,
      correlation = runif(1, -0.8, -0.5),  # Moderate to strong negative correlation
      direction = -direction,  # Opposite direction
      rand_number = rand_number
    )
    
  } else if (pattern_type == 2) {
    # Pattern 2: Strong correlation in responders, no correlation in non-responders
    
    # For responders: strong correlation (either positive or negative)
    corr_sign <- sample(c(-1, 1), 1)
    resp_data <- data_correlated_generation(
      n = length(responder_indices),
      mean_x = mean_val, mean_y = mean_val, 
      sd_x = sd_val, sd_y = sd_val,
      correlation = corr_sign * runif(1, 0.7, 0.9),  # Strong correlation
      direction = direction,
      rand_number = rand_number
    )
    
    # For non-responders: no correlation
    nonresp_data <- data_correlated_generation(
      n = length(nonresponder_indices),
      mean_x = mean_val, mean_y = mean_val, 
      sd_x = sd_val, sd_y = sd_val,
      correlation = runif(1, -0.2, 0.2),  # No substantial correlation
      direction = -direction,  # Opposite direction
      rand_number = rand_number
    )
    
  } else if (pattern_type == 3) {
    # Pattern 3: Opposite strong correlations (pos vs neg)
    
    # For responders: strong positive correlation
    resp_data <- data_correlated_generation(
      n = length(responder_indices),
      mean_x = mean_val, mean_y = mean_val, 
      sd_x = sd_val, sd_y = sd_val,
      correlation = runif(1, 0.6, 0.9),  # Strong positive correlation
      direction = direction,
      rand_number = rand_number
    )
    
    # For non-responders: strong negative correlation
    nonresp_data <- data_correlated_generation(
      n = length(nonresponder_indices),
      mean_x = mean_val, mean_y = mean_val, 
      sd_x = sd_val, sd_y = sd_val,
      correlation = runif(1, -0.9, -0.6),  # Strong negative correlation
      direction = -direction,  # Opposite direction
      rand_number = rand_number
    )
    
  } else {
    # Pattern 4: Different slopes in the correlation (strong vs weak same direction)
    corr_direction <- sample(c(-1, 1), 1)
    
    # For responders: strong correlation
    resp_data <- data_correlated_generation(
      n = length(responder_indices),
      mean_x = mean_val, mean_y = mean_val, 
      sd_x = sd_val, sd_y = sd_val,
      correlation = corr_direction * runif(1, 0.7, 0.9),  # Strong correlation
      direction = direction,
      rand_number = rand_number
    )
    
    # For non-responders: weak correlation in same direction
    nonresp_data <- data_correlated_generation(
      n = length(nonresponder_indices),
      mean_x = mean_val, mean_y = mean_val, 
      sd_x = sd_val, sd_y = sd_val,
      correlation = corr_direction * runif(1, 0.1, 0.3),  # Weak correlation
      direction = -direction,  # Opposite direction
      rand_number = rand_number
    )
  }
  
  # Replace values with correlated data
  result_matrix <- matrix_data
  result_matrix[entity1_idx, responder_indices] <- resp_data[, 1]
  result_matrix[entity2_idx, responder_indices] <- resp_data[, 2]
  result_matrix[entity1_idx, nonresponder_indices] <- nonresp_data[, 1]
  result_matrix[entity2_idx, nonresponder_indices] <- nonresp_data[, 2]
  
  return(result_matrix)
}

rnaseq_generation <- function(metadata) {
  # Use gene symbols from the metapathway
  gene_names <- all_genes
  
  # If we need more genes than in the metapathway, truncate to the number we need
  if (length(gene_names) > n_genes) {
    gene_names <- gene_names[1:n_genes]
  }
  
  # Base expression values (log2 normalized)
  rnaseq_matrix <- matrix(rnorm(length(gene_names) * n_patients, mean = 5, sd = 1), 
                          nrow = length(gene_names), ncol = n_patients)
  
  # Identify responders and non-responders
  responder_indices <- which(metadata$response == "Responder")
  nonresponder_indices <- which(metadata$response == "Non-responder")
  
  # Create gene pairs directly from the metapathway network
  gene_pairs <- list()
  pair_count <- 0
  
  # Create gene pairs from the metapathway
  for (i in 1:nrow(metapathway_gene_symbols)) {
    gene1 <- metapathway_gene_symbols$from[i]
    gene2 <- metapathway_gene_symbols$to[i]
    
    # Find indices of these genes in our gene_names vector
    gene1_idx <- which(gene_names == gene1)
    gene2_idx <- which(gene_names == gene2)
    
    # Only add if both genes are in our dataset
    if (length(gene1_idx) > 0 && length(gene2_idx) > 0) {
      pair_count <- pair_count + 1
      gene_pairs[[pair_count]] <- c(gene1_idx[1], gene2_idx[1])
    }
  }
  
  # # Report how many pairs were generated
  # cat(sprintf("Generated %d gene pairs in total\n", length(gene_pairs)))
  
  # For each gene pair, create differential correlation patterns 
  for (i in 1:length(gene_pairs)) {
    gene1_idx <- gene_pairs[[i]][1]
    gene2_idx <- gene_pairs[[i]][2]
    
    rnaseq_matrix <- correlation_differential_generation(
      responder_indices, nonresponder_indices,
      mean_val = 5, sd_val = 2,
      entity1_idx = gene1_idx, entity2_idx = gene2_idx,
      matrix_data = rnaseq_matrix
    )
  }
  
  # Ensure values stay in a realistic range (log2 normalized data typically ~0-20)
  rnaseq_matrix[rnaseq_matrix < 0] <- 0
  rnaseq_matrix[rnaseq_matrix > 23] <- 23
  
  # Set gene names
  rownames(rnaseq_matrix) <- gene_names
  colnames(rnaseq_matrix) <- metadata$patient_id
  
  return(rnaseq_matrix)
}

proteomics_generation <- function(metadata, rnaseq_matrix) {
  # Use gene names from RNA-seq as protein names where possible
  # This ensures consistency between omics layers
  gene_names <- rownames(rnaseq_matrix)
  protein_ids <- gene_names[1:min(length(gene_names), n_proteins)]
  
  # Base proteomics values (log2 normalized)
  proteomics_matrix <- matrix(rnorm(length(protein_ids) * n_patients, mean = 6,
                                    sd = 0.8), 
                              nrow = length(protein_ids), ncol = n_patients)
  
  # Identify responders and non-responders
  responder_indices <- which(metadata$response == "Responder")
  nonresponder_indices <- which(metadata$response == "Non-responder")
  
  n_diff_corr_pairs <- min(50, nrow(metapathway_gene_symbols))
  
  # Create protein pairs from the metapathway network
  protein_pairs <- list()
  pair_count <- 0

  for (i in 1:nrow(metapathway_gene_symbols)) {
    if (pair_count >= n_diff_corr_pairs) break
    
    protein1 <- metapathway_gene_symbols$from[i]
    protein2 <- metapathway_gene_symbols$to[i]
    
    protein1_idx <- which(protein_ids == protein1)
    protein2_idx <- which(protein_ids == protein2)
    
    # Only add if both proteins are in our dataset
    if (length(protein1_idx) > 0 && length(protein2_idx) > 0) {
      pair_count <- pair_count + 1
      protein_pairs[[pair_count]] <- c(protein1_idx[1], protein2_idx[1])
    }
  }
  
  # For each protein pair, create differential correlation patterns
  for (i in 1:length(protein_pairs)) {
    protein1_idx <- protein_pairs[[i]][1]
    protein2_idx <- protein_pairs[[i]][2]
    
    proteomics_matrix <- correlation_differential_generation(
      responder_indices, nonresponder_indices,
      mean_val = 6, sd_val = 1,
      entity1_idx = protein1_idx, entity2_idx = protein2_idx,
      matrix_data = proteomics_matrix
    )
  }
  
  # Ensure values stay in a realistic range for proteomics (typically 2-12 for log2 normalized data)
  proteomics_matrix[proteomics_matrix < 2] <- 2
  proteomics_matrix[proteomics_matrix > 12] <- 12
  
  # Set protein IDs and sample names
  rownames(proteomics_matrix) <- protein_ids
  colnames(proteomics_matrix) <- metadata$patient_id
  
  # Return just the matrix
  return(proteomics_matrix)
}

olink_generation <- function(metadata, rnaseq_matrix, proteomics_matrix) {
  # Use some protein names from proteomics data where possible to ensure consistency
  protein_ids <- rownames(proteomics_matrix)
  olink_ids <- protein_ids[1:min(length(protein_ids), n_olink)]
  
  # Base Olink values (NPX values typically range between 1-10)
  olink_matrix <- matrix(rnorm(length(olink_ids) * n_patients, mean = 5, sd = 1.2), 
                         nrow = length(olink_ids), ncol = n_patients)
  
  # Identify responders and non-responders
  responder_indices <- which(metadata$response == "Responder")
  nonresponder_indices <- which(metadata$response == "Non-responder")
  
  # Define number of differentially correlated protein pairs
  n_diff_corr_pairs <- min(15, nrow(metapathway_gene_symbols))
  
  # Create protein pairs from the metapathway network
  olink_pairs <- list()
  pair_count <- 0
  
  # Create Olink protein pairs from the metapathway
  for (i in 1:nrow(metapathway_gene_symbols)) {
    if (pair_count >= n_diff_corr_pairs) break
    
    protein1 <- metapathway_gene_symbols$from[i]
    protein2 <- metapathway_gene_symbols$to[i]
    
    # Find indices of these proteins in our olink_ids vector
    protein1_idx <- which(olink_ids == protein1)
    protein2_idx <- which(olink_ids == protein2)
    
    # Only add if both proteins are in our dataset
    if (length(protein1_idx) > 0 && length(protein2_idx) > 0) {
      pair_count <- pair_count + 1
      olink_pairs[[pair_count]] <- c(protein1_idx[1], protein2_idx[1])
    }
  }
  
  # For each Olink protein pair, create differential correlation patterns
  for (i in 1:length(olink_pairs)) {
    protein1_idx <- olink_pairs[[i]][1]
    protein2_idx <- olink_pairs[[i]][2]
    
    olink_matrix <- correlation_differential_generation(
      responder_indices, nonresponder_indices,
      mean_val = 5, sd_val = 1,
      entity1_idx = protein1_idx, entity2_idx = protein2_idx,
      matrix_data = olink_matrix
    )
  }
  
  # Ensure values stay in a realistic range for Olink (NPX values typically 1-10)
  olink_matrix[olink_matrix < 1] <- 1
  olink_matrix[olink_matrix > 10] <- 10
  
  # Set protein IDs and sample names
  rownames(olink_matrix) <- olink_ids
  colnames(olink_matrix) <- metadata$patient_id
  
  return(olink_matrix)
}

metadata <- metadata_generation()
rownames(metadata) <- metadata$patient_id
rnaseq <- rnaseq_generation(metadata)
proteomics <- proteomics_generation(metadata, rnaseq)
olink <- olink_generation(metadata, rnaseq, proteomics)

saveRDS(metadata, "data/synthetic_metadata.RDS")
saveRDS(rnaseq, "data/synthetic_rnaseqData.RDS")
saveRDS(proteomics, "data/synthetic_proteomicsData.RDS")
saveRDS(olink, "data/synthetic_OlinkData.RDS")
