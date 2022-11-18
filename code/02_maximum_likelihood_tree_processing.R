# metazoan-mixtures/code/02_maximum_likelihood_tree_processing.R
## This script processes maximum likelihood trees (so all trees from all datasets are consistent)
# Caitlin Cherryh, 2022

## This script:
# 1. Updates taxa names for maximum likelihood trees to be consistent across datasets



#### 1. Input parameters ####
## Specify parameters:
# iqtree_file_dir    <- Directory for IQ-Tree output (.log, .iqtree and .treefile files from IQ-Tree runs)
# renamed_tree_dir   <- Directory for results and plots
# repo_dir           <- Location of caitlinch/metazoan-mixtures github repository

iqtree_file_dir      <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/02_maximum_likelihood_trees/01_ml_tree_output_files/"
renamed_tree_dir     <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/02_maximum_likelihood_trees/02_renamed_trees/"
repo_dir        <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"



#### 2. Prepare functions, variables and packages ####
# Open packages
library(ape) # for `read.tree` and `write.tree`
library(TreeTools) # for `as.multiPhylo`

# Source functions and taxa lists
source(paste0(repo_dir, "code/func_naming.R"))

# Open the renaming csv
taxa_df <- read.csv(paste0(repo_dir, "Cherryh_MAST_metazoa_taxa_reconciliation.csv"), stringsAsFactors = FALSE)

#### 3. Update the taxa labels in each tree ####
# Extract the full list of trees
all_files <- paste0(iqtree_file_dir, list.files(iqtree_file_dir))
all_tree_files <- grep("\\.treefile", all_files, value = T)
# Rename all trees
all_trees_list <- lapply(all_tree_files, update.tree.taxa, naming_reconciliation_df = taxa_df, 
                    output.clade.names = TRUE, save.updated.tree = TRUE, output.directory = renamed_tree_dir)
all_trees <- as.multiPhylo(all_trees_list)


# Fixing Moroz taxa names
nosenko_trees_nonribo <- grep("nonribosomal",grep("Nosenko2013", all_tree_files, value = T), value = T)
treefile <- nosenko_trees_nonribo[1]

# Renaming function parameters
naming_reconciliation_df = taxa_df
output.clade.names = TRUE
save.updated.tree = TRUE
output.directory = renamed_tree_dir

# Find missing tips
tip_names <- read.tree(treefile)$tip.label
tip_names[which(!tip_names %in% tree_taxa_df$original_name)]
setdiff(t_names,tree_taxa_df$original_name)
# Check whether missing taxa are in mastmet_df
mastmet_file_path <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/Cherryh_MAST_metazoa_taxa_collation.csv"
mastmet_df <- read.csv(mastmet_file_path, stringsAsFactors = F)


