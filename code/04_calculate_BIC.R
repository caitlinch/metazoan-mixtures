## caitlinch/metazoan-mixtures/code/05_calculate_BIC.R
# This script calculates BIC values for each MAST analysis
# Caitlin Cherryh 2023


#### 1. Input parameters ####
## File paths
# repo_dir                <- Location of caitlinch/metazoan-mixtures github repository
# output_file_dir         <- Directory for output csvs
# hypothesis_tree_dir     <- Directory containing all constrained ML trees (i.e., the hypothesis trees)

## File paths
repo_dir                <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"
output_file_dir         <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/01_output_files/"
hypothesis_tree_dir     <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/04_hypothesis_trees/"



#### 2. Prepare variables, open packages and source functions ####
# Change the default number of digits (so we don't lose the decimal points in the BIC scores)
options(digits = 12)

# Open packages
library(ape)
library(phangorn)

# Source files containing functions
source(paste0(repo_dir, "code/func_data_processing.R"))

# List all files in output directory
all_output_files <- paste0(output_file_dir, list.files(output_file_dir))

# List all hypothesis trees
all_hypothesis_tree_paths <- paste0(hypothesis_tree_dir, list.files(hypothesis_tree_dir, recursive = T))



#### 1. Open csv files ####
# Open the BIC df
bic_df <- read.csv(grep("06_rank_BIC", all_output_files, value = T), stringsAsFactors = FALSE)
# Separate into tree and MAST rows
tree_bic_df <- bic_df[which(bic_df$num_trees == 1), ]
rownames(tree_bic_df) <- 1:nrow(tree_bic_df)
mast_bic_df <- bic_df[which(bic_df$num_trees > 1), ]
rownames(mast_bic_df) <- 1:nrow(mast_bic_df)



#### 2. Calculate the number of different branches for each MAST analysis ####
# check rows 185, 186
row_id <- 185
# Extract row
temp_row <- mast_bic_df[row_id, ]
# Extract hypothesis trees for this dataset/matrix/model class combination
temp_all_h_trees <- grep("\\.treefile", grep(temp_row$model_class, grep(temp_row$matrix, grep(temp_row$dataset, all_hypothesis_tree_paths, value = T), value = T), value = T), value = T)
# Extract tree files for this MAST run
if (temp_row$num_trees == 2){
  temp_h_trees <- c(grep("ML_H1", temp_all_h_trees, value = T), grep("ML_H2", temp_all_h_trees, value = T))
} else if (temp_row$num_trees == 3){
  temp_h_trees <- c(grep("ML_H1", temp_all_h_trees, value = T), grep("ML_H2", temp_all_h_trees, value = T), grep("ML_H3", temp_all_h_trees, value = T))
} else if (temp_row$num_trees == 5){
  temp_h_trees <- c(grep("ML_H1", temp_all_h_trees, value = T), grep("ML_H2", temp_all_h_trees, value = T), 
                    grep("ML_H3", temp_all_h_trees, value = T), grep("ML_H4", temp_all_h_trees, value = T), 
                    grep("ML_H5", temp_all_h_trees, value = T))
  # Compare splits in trees
}

trees_path <- temp_h_trees 

compare.splits.2.trees <- function(trees_path){
  ## Take 2 trees and calculate the number of different splits
  # Calculate the number of different splits using the RF distance
  t1 <- read.tree(trees_path[1])
  t2 <- read.tree(trees_path[2])
  num_different_splits <- RF.dist(t1, t2)
  return(num_different_splits)
}

compare.splits.3.trees <- function(trees_path){
  ## Take 3 trees and calculate the number of different splits
  # Open trees
  t1 <- read.tree(trees_path[1])
  t2 <- read.tree(trees_path[2])
  t3 <- read.tree(trees_path[2])
  # Calculate splits from trees and change into a matrix
  m1 <- as.data.frame(as.matrix(as.splits(t1)))
  m2 <- as.data.frame(as.matrix(as.splits(t2)))
  m3 <- as.data.frame(as.matrix(as.splits(t3)))
  # Rearrange columns to be in the same order
  col_order <- sort(colnames(m1))
  m1 <- m1[, col_order]
  m2 <- m2[, col_order]
  m3 <- m3[, col_order]
  # Combine all splits into a single dataframe
  m_all <- rbind(m1, m2, m3)
  # Identify splits that occur in all 3 trees
  duplicated_splits <- m_all |>
    group_by_all() |>
    filter(n() > 1) |>
    ungroup()
  unique_splits <- m_all |>
    group_by_all() |>
    filter(n() < 3) |>
    ungroup()
  # Check whether any splits in the unique_splits are replicates
  
  
  # Remove any splits in duplicated_splits that occur in unique_splits
  #   (first instance of a particular value not counted in duplicated)
  
  return(num_different_splits)
}



#### 3. Calculate BIC ####


