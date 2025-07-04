# Load necessary libraries
library(data.table)
library(matrixStats)
library(dplyr)
library(tibble)
library(stringr)
library(glue)
library(pheatmap)
library(grid)
library(JASPAR2022)
library(TFBSTools)
library("ComplexHeatmap")
library("circlize") #test

# Define variables
cancers <- c("breast", "bile_duct", "lung", "ovarian", "gastric", "pancreatic", "crc")
base_folder <- "result_TFs_diff_new_version/"
tfs_mat_dist <- matrix(0, nrow = 644, ncol = length(cancers))  # Fixed matrix column count
tfs_mat_corr <- matrix(0, nrow = 644, ncol = length(cancers))  # Fixed matrix column count
tfs_mat_dist_pre <- matrix(0, nrow = 644, ncol = length(cancers))  # Fixed matrix column count
tfs_mat_dist_post <- matrix(0, nrow = 644, ncol = length(cancers))  # Fixed matrix column count


jasp_names <- read.csv("JASPAR2022_ID_genes2.txt", sep = "\t", header = FALSE)
tf_rownames <-c()

# Loop through cancers and process files
for (i in seq_along(cancers)) {
  cancer <- cancers[i]
  print(cancer)
  folder <- glue("{base_folder}TF_{cancer}_diff/human/")
  flist <- list.files(folder, pattern = "_dist.tsv")
  names <-rep(NA, length(flist))
  species <- rep(NA, length(flist))
  
  for (j in seq_along(flist)) {
    f <- flist[j]
    indiv <- fread(paste0(folder, f), data.table = FALSE)
    
    # Ensure matrix dimensions are appropriate
    if (j <= nrow(tfs_mat_dist)) {
      tfs_mat_dist_pre[j, i] <- indiv$dist_distrib[970] - indiv$dist_distrib_2[970]
      tfs_mat_dist_post[j, i] <- indiv$dist_distrib[1030] - indiv$dist_distrib_2[1030]
      tfs_mat_dist[j, i] <- indiv$dist_distrib[1000] - indiv$dist_distrib_2[1000]
      tfs_mat_corr[j, i] <- cor(indiv$dist_distrib[900:1100],indiv$dist_distrib_2[900:1100])
    }
    
    # Extract gene names safely
    name <- gsub("\\..*", "", str_split(f, ".tsv", simplify = TRUE)[, 1]) 
    matched_name <- tryCatch(getMatrixByID(JASPAR2022, ID = name)@name, error=function(err) NA)
    names[j] <- ifelse(is.na(matched_name), "Unknown", matched_name)
    spe <- tryCatch(getMatrixByID(JASPAR2022, ID = name)@tags$species, error=function(err) NA)
    species[j] <- ifelse(is.na(spe), "Unknown", spe)
  }
}

rownames(tfs_mat_dist) <- names
colnames(tfs_mat_dist) <- cancers
spe_human <- which(species=="Homo sapiens")
tfs_mat_dist <- tfs_mat_dist[spe_human,]
tfs_mat_dist <- tfs_mat_dist[rownames(tfs_mat_dist != "Unkown"),]

rownames(tfs_mat_corr) <- names
colnames(tfs_mat_corr) <- cancers
spe_human <- which(species=="Homo sapiens")
tfs_mat_corr <- tfs_mat_corr[spe_human,]
tfs_mat_corr <- tfs_mat_corr[rownames(tfs_mat_corr != "Unkown"),]

scale_matrix <- function(mat) {
  max_abs_val <- max(abs(mat), na.rm = TRUE) # Find the maximum absolute value
  scaled_mat <- mat / max_abs_val # Scale while preserving distribution
  return(scaled_mat)
}

tfs_mat_dist_norm <- scale_matrix(tfs_mat_dist)

# Function to calculate quantiles for each column
calculate_quantiles <- function(data, probs = c(0.2)) {
  apply(data, 2, quantile, probs = probs)
}

# Calculate quantiles
med <- calculate_quantiles(abs(tfs_mat_dist_norm))

tfs_mat <- matrix(0, nrow = nrow(tfs_mat_dist_norm), ncol = ncol(tfs_mat_dist_norm))  # Fixed matrix column count
mat_features <- matrix(0, nrow = nrow(tfs_mat_dist_norm), ncol = ncol(tfs_mat_dist_norm)) 
vect_d_c <- 0
feat <- "Feature5"
# Loop through cancers and process files
for (i in (1:ncol(tfs_mat_dist_norm))) {
  for (j in (1:nrow(tfs_mat_dist_norm))) {
    
    dist <- tfs_mat_dist_norm[j,i]
    corre <- tfs_mat_corr[j,i]
    
    if (abs(dist) > abs(med[i]) && corre > 0.6 && dist < 0){
      vect_d_c <- abs(corre)
      feat <- "Feature1"
    }
    else if (abs(dist) > abs(med[i]) && corre < -0.6 && dist < 0){
      vect_d_c <- abs(corre)+10
      feat <- "Feature2"
    }
    else if (abs(dist) > abs(med[i]) && corre >0.6 && dist > 0){
      vect_d_c <- abs(corre) +20
      feat <- "Feature3"
    }
    else if (abs(dist) > abs(med[i]) && corre < -0.6 && dist > 0){
      vect_d_c <- abs(corre) +30
      feat <- "Feature4"
    }
    else {
      vect_d_c <- 0
      feat <- "Feature5"
    }
    tfs_mat[j,i] <- vect_d_c
    mat_features[j,i] <- feat
  }
}

rownames(tfs_mat) <- rownames(tfs_mat_dist_norm)
colnames(tfs_mat) <- cancers
rownames(mat_features) <- rownames(tfs_mat_dist_norm)
colnames(mat_features) <- cancers

col_fun_prop = colorRamp2(c(0, 0.5, 1,1.1, 10,10.5, 11, 11.1, 20, 20.5, 21, 21.1, 30,30.5,31, 31.1), 
                          c("white","darkgoldenrod1", "firebrick3","white", "white","lightblue1","#4169E1","white","white","palegreen", "darkolivegreen", "white","white","thistle2", "mediumpurple4","white"))

Heatmap(tfs_mat, col= col_fun_prop, width = unit(10, "cm"),show_column_dend = FALSE, show_row_dend = FALSE, show_row_names = F,border = TRUE)

