## caitlinch/metazoan-mixtures/code/05_calculate_BIC.R
# This script calculates BIC values for each MAST analysis
# Caitlin Cherryh 2023


#### 1. Input parameters ####
## File paths
# repo_dir              <- Location of caitlinch/metazoan-mixtures github repository
# output_file_dir       <- Directory for output csvs
# hypothesis_tree_dir   <- Directory containing all constrained ML trees (i.e., the hypothesis trees)

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
source(paste0(repo_dir, "code/func_BIC.R"))

# List all files in output directory
all_output_files          <- paste0(output_file_dir, list.files(output_file_dir))

# List all hypothesis trees
all_hypothesis_tree_paths <- paste0(hypothesis_tree_dir, list.files(hypothesis_tree_dir, recursive = T))
all_hypothesis_tree_paths <- grep("00_|2trees|5trees", all_hypothesis_tree_paths, value = T, invert = T)



#### 1. Open csv files ####
# Open the BIC df
in_df                 <- read.csv(grep("05_recalculate_BIC.csv", all_output_files, value = T), stringsAsFactors = FALSE)
new_col_order         <- grep("row_num|mast_nbranches|new_nfp|diff_nfp|new_BIC", names(in_df), value = T, invert = T)
in_df                 <- in_df[, new_col_order]
# Separate into tree and MAST rows
tree_bic_df           <- in_df[which(in_df$num_trees == 1), ]
rownames(tree_bic_df) <- 1:nrow(tree_bic_df)
mast_bic_df           <- in_df[which(in_df$num_trees > 1), ]
rownames(mast_bic_df) <- 1:nrow(mast_bic_df)



#### 2. Calculate the number of different branches for each MAST analysis ####
# Calculate number of +TR branches to consider per MAST analysis
tree_bic_df$num_unique_splits <- 0
tree_bic_df$num_shared_splits <- tree_bic_df$num_branches
split_output                  <- lapply(1:nrow(mast_bic_df), calculate.MAST.TR.branches, 
                                        MAST_output_df = mast_bic_df, 
                                        all_hypothesis_tree_paths = all_hypothesis_tree_paths)
mast_bic_df$num_unique_splits <- unlist(lapply(split_output, function(x){x[["num_unique_splits"]]}))
mast_bic_df$num_shared_splits <- unlist(lapply(split_output, function(x){x[["num_shared_splits"]]}))
# Subtract 1 from the number of shared splits in mast_bic_df
#     We need the number of splits in UNROOTED tree, whereas we calculated the number of splits in the ROOTED tree
mast_bic_df$num_shared_splits <- mast_bic_df$num_shared_splits - 1

# Collate dataframes
bic_df        <- rbind(tree_bic_df, mast_bic_df)
# Reorder dataframe
bic_df        <- bic_df[order(bic_df$dataset, bic_df$matrix, bic_df$model_class, bic_df$num_trees), ]
# Rearrange dataframe
bic_df        <- bic_df[ , c(1:19, 26:27, 20:25)]
names(bic_df) <- c(names(bic_df)[1:21], paste0("iqtree_", names(bic_df)[22:27]))



#### 3. Calculate BIC ####
# Calculate new number of free parameters
bic_df$new_num_free_params <- bic_df$model_np + bic_df$mixture_component_np + 
  bic_df$rhas_np + bic_df$state_freq_num_params + bic_df$num_MAST_weights + 
  bic_df$num_unique_splits + bic_df$num_shared_splits
# Compare number of free params 
# Note: Some MAST models were run using an old version of IQ-Tree where the number of free parameters were calculated using the +T model not the +TR model
#       therefore in some cases: diff_num_free_params>0
bic_df$diff_num_free_params <- abs(bic_df$new_num_free_params - bic_df$iqtree_num_free_params)
# Calculate BIC
bic_df$new_BIC              <- (-2*bic_df$iqtree_tree_LogL) + (bic_df$new_num_free_params * log(bic_df$nsites))
# Write file to csv
new_BIC_path                <- paste0(output_file_dir, "06_complete_BIC.csv")
write.csv(bic_df, file = new_BIC_path, row.names = FALSE)



#### 3. Create reduced table for results ####
# Add year column

# Factor model_class column for correct ordering (simple -> complex)

# Order by: year, dataset, matrix, model_class, num_trees


