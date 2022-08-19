# metazoan-mixtures/code/01_estimate_all_trees.R
## This script estimates maximum likelihood trees and maximum likelihood under constraint trees for empirical data sets
# Caitlin Cherryh, 2022

## This script:
# 1. Estimates ML trees for empirical data sets under different models
# 2. Estimates ML trees using a constraint tree for empirical data sets under different models
# 3. Applies the MAST (Mixtures Across Sites and Trees) model 



#### 1. Input parameters ####
## Specify parameters:
# alignment_dir     <- Directory containing alignments for all data sets
#                      Alignments have the naming convention dataset.matrix_name.sequence_type.fa
#                      E.g. Cherryh2022.all_taxa.aa.fa
# output_dir        <- Directory for IQ-Tree output (trees and tree mixtures)
# repo_dir          <- Location of caitlinch/metazoan-mixtures github repository

# iqtree2           <- Location of IQ-Tree2 stable release
# iqtree_tm         <- Location of IQ-Tree2 MAST release

location = "local"
if (location == "local"){
  alignment_dir <- "/Users/caitlin/Documents/PhD/Ch03_sponge_mixtures/01_alignments"
  output_dir <- "/Users/caitlin/Documents/PhD/Ch03_sponge_mixtures/02_output/"
  repo_dir <- "/Users/caitlin/Repositories/metazoan-mixtures/"
  
  iqtree2 <- "/Users/caitlin/Documents/PhD/Ch03_sponge_mixtures/iqtree-2.2.0-MacOSX/bin/iqtree2"
  iqtree2_tm <- "/Users/caitlin/Documents/PhD/Ch03_sponge_mixtures/iqtree-2.2.0.7.mix-MacOSX/bin/iqtree2"
  
} else if (location == "soma"){
  alignment_dir <- ""
  output_dir <- ""
  repo_dir <- ""
  
  iqtree2 <- ""
  iqtree2_tm <- ""
}



#### 2. Prepare functions, variables and packages ####
# Source functions
source(paste0(repo_dir, "code/func_constraint_trees.R"))
source(paste0(repo_dir, "code/func_data_processing.R"))

# Source information about datasets
source(paste0(repo_dir, "code/data_dataset_info.R"))

# Sort and categorise models of sequence evolution 
all_models <- sort(unique(unlist(lapply(all_models, sort.model.chunks))))
# Identify unique matrices within the models
model_components <- sort(unique(unlist(strsplit(all_models, "\\+"))))
# Remove components that do not correspond to a matrix or model (e.g. remove "I", "F", "G4")
model_components <- model_components[!model_components %in% c("F", "FO", "G", "G4", "I", "R", "R4")]
# PMSF is more complex than other models - requires separate function
# For now, remove PMSF from list 
model_components <- model_components[!model_components == "PMSF"]
# Note: partitioning schemes currently not possible in mixture of trees implementation



#### 3. Process each dataset for each model of sequence evolution ####
# Extract the list of alignments
all_alignments <- list.files(alignment_dir, full.names = T)

# For each alignment:
for (a in all_alignments){
  # Extract details about alignment from file name
  # Identify dataset (original paper citation)
  a_dataset <- strsplit(basename(a), "\\.")[[1]][1]
  # Identify matrix (which alignment file from original reference)
  a_matrix_id <- strsplit(basename(a), "\\.")[[1]][2]
}


# - Take a single alignment:
#   - Take a single model of sequence evolution
#     - Use ModelFinder to determine the best model within that model category
#     - Estimate an ML tree with that model of sequuence evolution
#     - Estimate the constraint trees
#     - Apply the mixture of trees method using +TR
#     - Apply the mixture of trees method using +T