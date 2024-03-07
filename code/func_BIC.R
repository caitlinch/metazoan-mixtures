## caitlinch/metazoan-mixtures/code/func_BIC.R
# Functions for manipulating, processing, and preparing datasets
# Caitlin Cherryh 2023


#### Packages ####
library(ape) # read.tree
library(phangorn) # as.splits
library(dplyr) # manipulating dataframes



#### Calculating number of different splits ####
calculate.MAST.TR.branches <- function(row_id, MAST_output_df, all_hypothesis_tree_paths){
  ## Calculate the number of splits for the MAST +TR model
  ##      i.e., the number of splits that occur in one or more tree that are NOT present in all trees
  
  # Extract row
  temp_row <- MAST_output_df[row_id, ]
  # Extract hypothesis trees for this dataset/matrix/model class combination
  temp_all_h_trees <- grep("\\.treefile", grep(temp_row$model_class, grep(temp_row$matrix, grep(temp_row$dataset, all_hypothesis_tree_paths, value = T), value = T), value = T), value = T)
  # Extract correct number of trees for this analysis
  if (temp_row$num_trees == 2){
    temp_h_trees <- c(grep("ML_H1", temp_all_h_trees, value = T), grep("ML_H2", temp_all_h_trees, value = T))
  } else if (temp_row$num_trees == 3){
    temp_h_trees <- c(grep("ML_H1", temp_all_h_trees, value = T), grep("ML_H2", temp_all_h_trees, value = T),
                      grep("ML_H3", temp_all_h_trees, value = T))
  } else if (temp_row$num_trees == 5){
    temp_h_trees <- c(grep("ML_H1", temp_all_h_trees, value = T), grep("ML_H2", temp_all_h_trees, value = T),
                      grep("ML_H3", temp_all_h_trees, value = T), grep("ML_H4", temp_all_h_trees, value = T),
                      grep("ML_H5", temp_all_h_trees, value = T))
  }
  # Compare splits in trees
  num_splits <- compare.splits.n.trees(trees_path = temp_h_trees)
  # Return output
  return(num_splits)
}



compare.splits.n.trees <- function(trees_path){
  ## Take n trees and calculate the number of shared/different splits
  
  ## Identify shared and unique splits
  # Find number of trees
  n_trees <- length(trees_path)
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
    filter(n() == n_trees) |>
    ungroup()
  unique_splits <- m_all |>
    group_by_all() |>
    filter(n() < n_trees) |>
    ungroup()
  
  ## Count number of unique splits
  # Check whether any splits in the unique_splits are replicates
  unique_unique_splits <- which(duplicated(unique_splits) == FALSE)
  duplicated_unique_splits <- which(duplicated(unique_splits) == TRUE)
  # Count number of unique splits
  # i.e., those splits where:
  #       - the split is present in less than five trees
  #       - the split is the first occurrence of that split - duplicates are not counted
  num_different_splits <- length(unique_unique_splits)
  
  ## Count number of shared splits
  # Identify duplicated splits in the all_tree_splits df
  shared_splits <- which(duplicated(all_tree_splits) == FALSE)
  # Count number of shared splits
  num_shared_splits <- length(shared_splits)
  
  ## Output
  op <- c(num_different_splits, num_shared_splits)
  names(op) <- c("num_unique_splits", "num_shared_splits")
  # Return the number of different splits
  return(op)
}



compare.splits.2.trees <- function(trees_path){
  ## Take 2 trees and calculate the number of different splits
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
    filter(n() == 2) |>
    ungroup()
  unique_splits <- m_all |>
    group_by_all() |>
    filter(n() < 2) |>
    ungroup()
  # Check whether any splits in the unique_splits are replicates
  unique_unique_splits <- which(duplicated(unique_splits) == FALSE)
  duplicated_unique_splits <- which(duplicated(unique_splits) == TRUE)
  # Count number of unique splits
  # i.e., those splits where:
  #       - the split is present in less than five trees
  #       - the split is the first occurrence of that split - duplicates are not counted
  num_different_splits <- length(unique_unique_splits)
  num_identical_splits <- nrow(all_tree_splits)
  # Return the number of different splits
  return(num_different_splits)
}



compare.splits.3.trees <- function(trees_path){
  ## Take 2 trees and calculate the number of different splits
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
  num_identical_splits <- nrow(all_tree_splits)
  # Return the number of different splits
  return(c(num_different_splits, num_identical_splits))
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
  num_identical_splits <- nrow(all_tree_splits)
  # Return the number of different splits
  return(c(num_different_splits, num_identical_splits))
}



#### Find best BIC for each analysis ####
find.best.bic <- function(row_id, params_df, bic_df){
  ## Find best BIC for each analysis
  
  # Extract parameters for this row
  temp_row <- params_df[row_id, ]
  # Extract rows from bic_df with these parameters
  temp_bic_df <- bic_df[which(bic_df$dataset == temp_row$dataset &
                                bic_df$matrix_name == temp_row$matrix_name &
                                bic_df$model_class == temp_row$model_class), ]
  # Find minimum BIC score
  min_BIC <- min(temp_bic_df$new_BIC)
  # Assemble output
  bic_output <- rep(min_BIC, nrow(temp_bic_df))
  # Return output
  return(bic_output)
}



