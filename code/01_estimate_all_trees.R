# metazoan-mixtures/code/01_estimate_all_trees.R
## This script estimates maximum likelihood trees and maximum likelihood under constraint trees for empirical data sets
# Caitlin Cherryh, 2022

## This script:
# 1. Estimates ML trees for empirical data sets under different models
# 2. Estimates ML trees using a constraint tree for empirical data sets under different models
# 3. Applies the MAST (Mixtures Across Sites and Trees) model 



#### 1. Input parameters ####
## Specify parameters:
# alignment_dir       <- Directory containing alignments for all data sets
#                        Alignments have the naming convention dataset.matrix_name.sequence_type.fa
#                        E.g. Cherryh2022.all_taxa.aa.fa
# output_dir          <- Directory for IQ-Tree output (trees and tree mixtures)
# repo_dir            <- Location of caitlinch/metazoan-mixtures github repository

# iqtree2             <- Location of IQ-Tree2 stable release
# iqtree_tm           <- Location of IQ-Tree2 MAST release

# iqtree_num_threads  <- Number of parallel threads for IQ-Tree to use. Can be a set number (e.g. 2) or "AUTO"

location = "local"
if (location == "local"){
  alignment_dir <- "/Users/caitlin/Documents/PhD/Ch03_sponge_mixtures/01_alignments"
  output_dir <- "/Users/caitlin/Documents/PhD/Ch03_sponge_mixtures/02_output/"
  repo_dir <- "/Users/caitlin/Repositories/metazoan-mixtures/"
  
  iqtree2 <- "/Users/caitlin/Documents/PhD/Ch03_sponge_mixtures/iqtree-2.2.0-MacOSX/bin/iqtree2"
  iqtree2_tm <- "/Users/caitlin/Documents/PhD/Ch03_sponge_mixtures/iqtree-2.2.0.7.mix-MacOSX/bin/iqtree2"
  
  iqtree_num_threads <- "AUTO"
  
} else if (location == "soma"){
  alignment_dir <- ""
  output_dir <- ""
  repo_dir <- ""
  
  iqtree2 <- ""
  iqtree2_tm <- ""
  
  iqtree_num_threads <- "AUTO"
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
# Extract the list of all files from the folder containing alignments/models/partition files
all_files <- list.files(alignment_dir, full.names = T)
# Extract the list of alignments (i.e. files that contain the word "alignment")
all_alignments <- grep("alignment", all_files, value = T)

# For each alignment:
for (a in all_alignments){
  # Extract details about alignment from file name
  # Identify dataset (original paper citation)
  a_dataset <- strsplit(basename(a), "\\.")[[1]][1]
  # Identify matrix (which alignment file from original reference)
  a_matrix_id <- strsplit(basename(a), "\\.")[[1]][2]
  
  # Create output folders for this dataset
  # Create a new folder in the output directory for this dataset
  a_op_dir <- paste0(output_dir, a_dataset, "/")
  if (dir.exists(a_op_dir) == FALSE){dir.create(a_op_dir)}
  # Create a new folder for storing maximum likelihood trees for this dataset
  a_ml_op_dir <- paste0(a_op_dir, "ML_trees", "/")
  if (dir.exists(a_ml_op_dir) == FALSE){dir.create(a_ml_op_dir)}
  # Create a new folder for storing constraint trees for this dataset
  a_c_op_dir <- paste0(a_op_dir, "constraint_trees", "/")
  if (dir.exists(a_c_op_dir) == FALSE){dir.create(a_c_op_dir)}
  # Create a new folder in the constraint folder for output for this dataset
  a_tm_op_dir <- paste0(a_op_dir, "tree_mixtures", "/")
  if (dir.exists(a_tm_op_dir) == FALSE){dir.create(a_tm_op_dir)}
  
  # For each of the model components:
  for (m in model_components){
    # Change directory to the maximum likelihood tree output directory for this dataset
    # This ensures IQ-Tree output files will be stored in the correct directory
    setwd(a_ml_op_dir)
    
    # Set a prefix for the ML tree for this combination of dataset, matrix, and model
    a_m_prefix <- paste0(a_dataset, ".", a_matrix_id, ".", m)
    
    # Estimate a maximum likelihood tree with the best model of sequence evolution containing that model component
    estimate.ml.iqtree(iqtree2, a, model = NA, mset = m, partition_file = NA, 
                       prefix = paste0(a_m_prefix, ".ML"), number_parallel_threads = "AUTO", number_of_bootstraps = NA,
                       redo = FALSE, safe = FALSE)
    
    # Extract information about this dataset
    a_info <- all_datasets[[a_dataset]]
    
    # Change directory to the constraint tree output directory for this dataset
    # This ensures IQ-Tree output files will be stored in the correct directory
    setwd(a_c_op_dir)
    
    # Create constraint trees
    constraint_df <- create.constraint.trees(dataset = a_dataset, tree_id = a_m_prefix, 
                                             dataset_constraint_tree_dir = a_c_op_dir, 
                                             model = m, model_id = m, outgroup_taxa = a_info$Outgroup,
                                             ctenophora_taxa = a_info$Ctenophora, porifera_taxa = a_info$Porifera, 
                                             sponges_1_taxa = as.character(unlist(a_info[c(a_info$Sponges_1)])), 
                                             sponges_2_taxa = as.character(unlist(a_info[c(a_info$Sponges_2)])), 
                                             placozoa_taxa = a_info$Placozoa, cnidaria_taxa = a_info$Cnidaria, 
                                             bilateria_taxa = a_info$Bilateria, alignment_file = a, 
                                             partitioned_check = FALSE, partition_file = NA, 
                                             iqtree_path = iqtree2, number_parallel_threads = iqtree_num_threads)
    
    # Estimate hypothesis trees for each of the constraint trees
    lapply(1:nrow(constraint_df), run.one.constraint.tree, constraint_df)
    
    # Combine hypothesis trees into one file and save
    hyp_tree_files <- combine.hypothesis.trees(tree_id = a_m_prefix, constraint_tree_directory = a_c_op_dir, 
                                               outgroup_taxa = a_info$Outgroup)
    
  } # end for (m in model_components)
} # end for (a in all_alignments)


# - Take a single alignment:
#   - Take a single model of sequence evolution
#     - Use ModelFinder to determine the best model within that model category
#     - Estimate an ML tree with that model of sequuence evolution
#     - Estimate the constraint trees
#     - Apply the mixture of trees method using +TR
#     - Apply the mixture of trees method using +T