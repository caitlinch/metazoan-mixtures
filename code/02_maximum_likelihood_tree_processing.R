# metazoan-mixtures/code/02_maximum_likelihood_tree_processing.R
## This script processes maximum likelihood trees (so all trees from all datasets are consistent)
# Caitlin Cherryh, 2022

## This script:
# 1. Plots figures for the introduction and methods section of the manuscript



#### 1. Input parameters ####
## Specify parameters:
# output_dir    <- Directory for IQ-Tree output (.log, .iqtree and .treefile files from IQ-Tree runs)
# results_dir   <- Directory for results and plots
# repo_dir      <- Location of caitlinch/metazoan-mixtures github repository

output_dir      <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/02_maximum_likelihood_trees/01_ml_tree_output_files/"
results_dir     <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/02_maximum_likelihood_trees/02_renamed_trees/"
repo_dir        <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"



#### 2. Prepare functions, variables and packages ####
# Open packages
library(ape)

# Source functions and taxa lists
source(paste0(repo_dir, "code/func_naming.R"))

# Open the renaming csv
naming_reconciliation_df <- read.csv(paste0(repo_dir, "MAST_metazoa_taxa_reconciliation.csv"), stringsAsFactors = FALSE)

#### 3. Update the taxa labels in each tree ####
# Extract the full list of trees
all_files <- paste0(output_dir, list.files(output_dir))
all_tree_files <- grep("\\.treefile", all_files, value = T)

treefile <- all_tree_files[1]
output.directory = results_dir
save.updated.tree = TRUE

