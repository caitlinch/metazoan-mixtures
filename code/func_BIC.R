## caitlinch/metazoan-mixtures/code/func_BIC.R
# Functions for manipulating, processing, and preparing datasets
# Caitlin Cherryh 2023


#### Packages ####
library(ape) # read.tree
library(phangorn) # as.splits
library(dplyr) # manipulating dataframes


#### Calculating number of different splits ####
calculate.MAST.TR.branches <- function(row_id, MAST_output_df){
  ## Calculate the number of splits for the MAST +TR model
  ##      i.e., the number of splits that occur in one or more tree that are NOT present in all trees
  
  # Extract row
  temp_row <- mast_bic_df[row_id, ]
  # Extract hypothesis trees for this dataset/matrix/model class combination
  temp_all_h_trees <- grep("\\.treefile", grep(temp_row$model_class, grep(temp_row$matrix, grep(temp_row$dataset, all_hypothesis_tree_paths, value = T), value = T), value = T), value = T)
  # Extract tree files for this MAST run and compare splits in trees
  if (temp_row$num_trees == 2){
    temp_h_trees <- c(grep("ML_H1", temp_all_h_trees, value = T), grep("ML_H2", temp_all_h_trees, value = T))
    num_diff_splits <- compare.splits.2.trees(trees_path = temp_h_trees)
  } else if (temp_row$num_trees == 3){
    temp_h_trees <- c(grep("ML_H1", temp_all_h_trees, value = T), grep("ML_H2", temp_all_h_trees, value = T), grep("ML_H3", temp_all_h_trees, value = T))
    num_diff_splits <- compare.splits.3.trees(trees_path = temp_h_trees)
  } else if (temp_row$num_trees == 5){
    temp_h_trees <- c(grep("ML_H1", temp_all_h_trees, value = T), grep("ML_H2", temp_all_h_trees, value = T), 
                      grep("ML_H3", temp_all_h_trees, value = T), grep("ML_H4", temp_all_h_trees, value = T), 
                      grep("ML_H5", temp_all_h_trees, value = T))
    num_diff_splits <- compare.splits.5.trees(trees_path = temp_h_trees)
  }
  # Return output
  return(num_diff_splits)
}



compare.splits.2.trees <- function(trees_path){
  ## Take 2 trees and calculate the number of different splits
  # Calculate the number of different splits using the RF distance
  t1 <- read.tree(trees_path[1])
  t2 <- read.tree(trees_path[2])
  num_different_splits <- RF.dist(t1, t2)
  return(num_different_splits)
}



compare.splits.3.trees <- function(trees_path){
  # Open trees
  t_list <- lapply(trees_path, read.tree)
  # Calculate splits from trees and change into a matrix
  m_list <- lapply(t_list, function(x){as.data.frame(as.matrix(as.splits(x)))})
  # Rearrange columns to be in the same order
  col_order <- sort(colnames(m_list[[1]]))
  m_ordered_list <- lapply(m_list, function(x){x[, col_order]})
  # Combine all splits into a single dataframe
  m_all <- bind_rows(m_ordered_list)
  # Identify splits that occur in all 5 trees
  all_tree_splits <- m_all |>
    group_by_all() |>
    filter(n() == 3) |>
    ungroup()
  unique_splits <- m_all |>
    group_by_all() |>
    filter(n() < 3) |>
    ungroup()
  # Check whether any splits in the unique_splits are replicates
  unique_unique_splits <- which(duplicated(unique_splits) == FALSE)
  duplicated_unique_splits <- which(duplicated(unique_splits) == TRUE)
  # Count number of unique splits
  # i.e., those splits where:
  #       - the split is present in less than five trees
  #       - the split is the first occurrence of that split - duplicates are not counted
  num_different_splits <- length(unique_unique_splits)
  # Return the number of different splits
  return(num_different_splits)
}



compare.splits.5.trees <- function(trees_path){
  ## Take 5 trees and calculate the number of different splits
  # Open trees
  t_list <- lapply(trees_path, read.tree)
  # Calculate splits from trees and change into a matrix
  m_list <- lapply(t_list, function(x){as.data.frame(as.matrix(as.splits(x)))})
  # Rearrange columns to be in the same order
  col_order <- sort(colnames(m_list[[1]]))
  m_ordered_list <- lapply(m_list, function(x){x[, col_order]})
  # Combine all splits into a single dataframe
  m_all <- bind_rows(m_ordered_list)
  # Identify splits that occur in all 5 trees
  all_tree_splits <- m_all |>
    group_by_all() |>
    filter(n() == 5) |>
    ungroup()
  unique_splits <- m_all |>
    group_by_all() |>
    filter(n() < 5) |>
    ungroup()
  # Check whether any splits in the unique_splits are replicates
  unique_unique_splits <- which(duplicated(unique_splits) == FALSE)
  duplicated_unique_splits <- which(duplicated(unique_splits) == TRUE)
  # Count number of unique splits
  # i.e., those splits where:
  #       - the split is present in less than five trees
  #       - the split is the first occurrence of that split - duplicates are not counted
  num_different_splits <- length(unique_unique_splits)
  # Return the number of different splits
  return(num_different_splits)
}



