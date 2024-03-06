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
library(readxl)

# Source files containing functions
source(paste0(repo_dir, "code/func_data_processing.R"))

# List all files in output directory
all_output_files <- paste0(output_file_dir, list.files(output_file_dir))
# # Remove any files with all 5 trees - only want to look at the output for the 2 tree model
# all_output_files <- grep("5trees", all_output_files, value = TRUE, invert = TRUE)


#### 10. Determine best BIC for each analysis ####
