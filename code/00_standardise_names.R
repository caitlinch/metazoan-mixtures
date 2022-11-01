# metazoan-mixtures/code/00_standardise_names.R
## This script standardises the names for all species for all datasets
# Caitlin Cherryh, 2022

# For this script, you will need the naming csv from Li et. al. (2021)
#     Available from the data repository: https://figshare.com/articles/dataset/Rooting_the_animal_tree_of_life/13122881
#     Download the reconciliation.tar.xz and identify the location of the "reconciliation/taxonomy_info/taxon_table.tsv" file



#### 1. Input parameters ####
## Specify parameters:
# output_dir          <- Directory for IQ-Tree output (trees and tree mixtures)
# repo_dir            <- Location of caitlinch/metazoan-mixtures github repository
# taxon_table_path    <- Location of the reconciliation/taxonomy_info/taxon_table.tsv from Li et. al. (2021)

location = "local"
if (location == "local"){
  output_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/"
  repo_dir <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"
  taxon_table_path <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/00_Li2021_supp/reconciliation_keep/taxonomy_info/taxon_table.tsv"
} 



#### 2. Open packages and source functions ####
# Source functions
source(paste0(repo_dir, "code/func_naming.R"))
source(paste0(repo_dir, "code/data_dataset_info.R"))



#### 3. Prepare csv ####



