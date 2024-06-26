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
library(stringr)

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
# Each tree has 0 "unique" splits and number of "shared" splits equal to the number of branches in the unrooted tree
#     Here, "shared" means present in all 1 tree
tree_bic_df$num_unique_splits <- 0
tree_bic_df$num_shared_splits <- tree_bic_df$num_branches
# Calculate number of +TR branches to consider per MAST analysis
split_output                  <- lapply(1:nrow(mast_bic_df), calculate.MAST.TR.branches, 
                                        MAST_output_df = mast_bic_df, 
                                        all_hypothesis_tree_paths = all_hypothesis_tree_paths)
# Add number of distinct branches in each MAST analysis as a column
mast_bic_df$num_unique_splits <- unlist(lapply(split_output, function(x){x[["num_unique_splits"]]}))
# Subtract 1 from the number of shared splits in mast_bic_df
#     We need the number of splits in UNROOTED tree, whereas we calculated the number of splits in the ROOTED tree
mast_bic_df$num_shared_splits <- unlist(lapply(split_output, function(x){x[["num_shared_splits"]]})) - 1

# Collate dataframes
bic_df  <- rbind(tree_bic_df, mast_bic_df)
# Reorder dataframe
bic_df  <- bic_df[order(bic_df$dataset, bic_df$matrix, bic_df$model_class, bic_df$num_trees), ]
# Rearrange dataframe
bic_df            <- bic_df[ , c(1:19, 26:27, 20:25)]
names(bic_df)     <- c(names(bic_df)[1:21], paste0("iqtree_", names(bic_df)[22:27]))
rownames(bic_df)  <- 1:nrow(bic_df)

# Create a new dataframe with the unique parameters for each analysis: dataset, matrix, and model
params_df           <- bic_df[ c("dataset", "matrix", "model_class")]
params_df           <- params_df[which(duplicated(params_df) == FALSE), ]
rownames(params_df) <- 1:nrow(params_df)



#### 3. Calculate BIC ####
# Calculate new number of free parameters
bic_df$new_num_free_params <- bic_df$model_np + bic_df$mixture_component_np + 
  bic_df$rhas_np + bic_df$state_freq_num_params + bic_df$num_MAST_weights + 
  bic_df$num_unique_splits + bic_df$num_shared_splits
# Compare number of free params 
# Note: Some MAST models were run using an old version of IQ-Tree where the number of free parameters were calculated using the +T model not the +TR model
#       therefore in some cases: diff_num_free_params>0
bic_df$diff_num_free_params <- abs(bic_df$new_num_free_params - bic_df$iqtree_num_free_params)
# Calculate missing branches
bic_df$excluded_duplicate_branches <- (bic_df$num_trees * bic_df$num_branches) - (bic_df$num_unique_splits + bic_df$num_shared_splits)
# Calculate BIC
bic_df$new_BIC  <- (-2*bic_df$iqtree_tree_LogL) + (bic_df$new_num_free_params * log(bic_df$nsites))
# Find minimum BIC for each analysis
bic_df$best_BIC <- unlist(lapply(1:nrow(params_df), find.best.bic, 
                                 params_df = params_df, 
                                 bic_df = bic_df))
# Calculate delta BIC for each analysis
bic_df$delta_BIC <- bic_df$new_BIC - bic_df$best_BIC

## Format and save dataframe
# Round newly calculated columns to 3 dp
bic_df$new_BIC    <- round(bic_df$new_BIC, digits = 3)
bic_df$best_BIC   <- round(bic_df$best_BIC, digits = 3)
bic_df$delta_BIC  <- round(bic_df$delta_BIC, digits = 3)
# Add year column
bic_df$year <- as.numeric(str_extract(bic_df$dataset, "(\\d)+"))
# Write file to csv
new_BIC_path    <- paste0(output_file_dir, "MS_complete_BIC.csv")
write.csv(bic_df, file = new_BIC_path, row.names = FALSE)



#### 3. Create reduced table for results ####
# Create new df for formatted results
pretty_bic_df <- bic_df[ , c("year", "dataset", "matrix", "model_class", "model", 
                             "num_trees", "tree_topology", "iqtree_tree_LogL", 
                             "new_num_free_params", "new_BIC", "delta_BIC")]
# Factor model_class column for correct ordering (simple -> complex)
# Levels: Single (Q), Mixture (Mixture), Posterior Mean Site Frequency (PMSF), Empirical Profile Mixture Model (PM)
pretty_bic_df$model_class <- factor(pretty_bic_df$model_class,
                                    levels = c("Single", "Other", "PMSF", "CXX"),
                                    labels = c("Q", "Mixture", "PMSF", "PM"),
                                    ordered = TRUE)
# Order by: year, dataset, matrix, model_class, num_trees
pretty_bic_df <- pretty_bic_df[order(pretty_bic_df$year, pretty_bic_df$dataset, pretty_bic_df$matrix,
                                     pretty_bic_df$model_class, pretty_bic_df$num_trees) , ]
pretty_bic_df <- pretty_bic_df[ , 2:ncol(pretty_bic_df)]
# Save as output
pretty_BIC_path    <- paste0(output_file_dir, "summary_all_BIC.csv")
write.csv(pretty_bic_df, file = pretty_BIC_path, row.names = FALSE)



#### 4. Output only best BIC for each analysis ####
# Create new df for formatted results
best_bic_df <- pretty_bic_df[which(pretty_bic_df$delta_BIC == 0), ]
# Add year column
best_bic_df$year <- as.numeric(str_extract(best_bic_df$dataset, "(\\d)+"))
# Order by: year, dataset, matrix, model_class, num_trees
best_bic_df <- best_bic_df[order(best_bic_df$year, best_bic_df$dataset, best_bic_df$matrix,
                                 best_bic_df$model_class, best_bic_df$num_trees) , ]
# Remove columns
best_bic_df <- best_bic_df[ , c(1:(ncol(best_bic_df) - 2))]
# Save as output
best_BIC_path    <- paste0(output_file_dir, "summary_best_BIC.csv")
write.csv(best_bic_df, file = best_BIC_path, row.names = FALSE)





