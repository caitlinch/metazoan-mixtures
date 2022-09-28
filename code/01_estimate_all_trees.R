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
  alignment_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/01_Data_all"
  output_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/"
  repo_dir <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"
  
  iqtree2 <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_Software_IQ-Tree/iqtree-2.2.0-MacOSX/bin/iqtree2"
  iqtree2_tm <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_Software_IQ-Tree/iqtree-2.2.0.7.mix-MacOSX/bin/iqtree2"
  
  iqtree_num_threads <- "AUTO"
  
} else if (location == "server"){
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
# For now, remove ModelFinder from the list (comparing existing models - ModelFinder will use best model)
model_components <- model_components[!model_components == "ModelFinder"]
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
                       redo = FALSE, safe = FALSE, run.iqtree = FALSE)
    
    # Extract the .iqtree file for the prefix
    all_ml_op_files <- paste0(a_ml_op_dir, list.files(a_ml_op_dir))
    all_iqtree_files <- grep("\\.iqtree", all_ml_op_files, value = TRUE)
    prefix_iqtree_file <- grep(a_m_prefix, all_iqtree_files, value = TRUE)
    
    # Extract best model and feed into analysis (estimating hypothesis trees; applying mixture of trees method)
    best_model <- extract.best.model(prefix_iqtree_file)
    
    # Extract information about model parameters
    model_parameters <- extract.model.details(prefix_iqtree_file)
    # Save model parameters as csv
    write.csv(model_parameters$parameters, file = paste0(a_ml_op_dir, a_m_prefix, "_model_parameters.csv"), row.names = FALSE)
    # Save gamma categories as csv (if present in .iqtree file)
    if (is.na(model_parameters$gamma_categories) == FALSE){
      write.csv(model_parameters$gamma_categories, file = paste0(a_ml_op_dir, a_m_prefix, "_model_gamma.categories.csv"), row.names = FALSE)
    }
    # Save state frequencies as csv (if present in .iqtree file)
    if (is.na(model_parameters$frequency) == FALSE){
      write.csv(model_parameters$frequency, file = paste0(a_ml_op_dir, a_m_prefix, "_model_state.frequencies.csv"), row.names = FALSE)
    }
    
    # Extract information about this dataset
    a_info <- all_datasets[[a_dataset]]
    
    # Change directory to the constraint tree output directory for this dataset
    # This ensures IQ-Tree output files will be stored in the correct directory
    setwd(a_c_op_dir)
    
    # Create constraint trees, including best model in the dataframe (so it will be used to estimate hypothesis trees)
    constraint_df <- create.constraint.trees(dataset = a_dataset, tree_id = a_m_prefix, 
                                             dataset_constraint_tree_dir = a_c_op_dir, 
                                             model = best_model, model_id = m, outgroup_taxa = a_info$Outgroup,
                                             ctenophora_taxa = a_info$Ctenophora, porifera_taxa = a_info$Porifera, 
                                             sponges_1_taxa = as.character(unlist(a_info[c(a_info$Sponges_1)])), 
                                             sponges_2_taxa = as.character(unlist(a_info[c(a_info$Sponges_2)])), 
                                             placozoa_taxa = a_info$Placozoa, cnidaria_taxa = a_info$Cnidaria, 
                                             bilateria_taxa = a_info$Bilateria, alignment_file = a, 
                                             partitioned_check = FALSE, partition_file = NA, 
                                             iqtree_path = iqtree2, number_parallel_threads = iqtree_num_threads,
                                             run.iqtree = FALSE)
    
    # Estimate hypothesis trees for each of the constraint trees
    lapply(1:nrow(constraint_df), run.one.constraint.tree, constraint_df)
    
    # Combine hypothesis trees into one file and save
    hyp_tree_files <- combine.hypothesis.trees(tree_id = a_m_prefix, constraint_tree_directory = a_c_op_dir, 
                                               outgroup_taxa = a_info$Outgroup)
    # Get name for rooted hypothesis trees
    rooted_hyp_trees <- hyp_tree_files[grep("rooted", names(hyp_tree_files))]
    
    # Change directory to the tree mixtures output directory for this dataset
    # This ensures IQ-Tree output files will be stored in the correct directory
    setwd(a_tm_op_dir)
    
    # Apply mixture of trees method with best model from maximum likelihood tree estimation
    # Run with +TR option (same branch lengths) 
    tree_branch_model <- "TR"
    run.tree.mixture.model(alignment_file = a, hypothesis_tree_file = rooted_hyp_trees, 
                           partition_file = NA, use.partition = FALSE, prefix = paste0(a_m_prefix,".", tree_branch_model),
                           model = best_model, iqtree2_tree_mixtures_implementation = iqtree2_tm, 
                           tree_branch_option = tree_branch_model, number_parallel_threads = iqtree_num_threads,
                           run.iqtree = FALSE)
    
    # Identify iqtree files from tree mixture run
    all_files <- paste0(a_tm_op_dir, list.files(a_tm_op_dir))
    all_iqtree_files <- grep("\\.iqtree", all_files, value = TRUE)
    tree_mixture_tr_iqfile <- grep(paste0(a_m_prefix,".TR"), all_iqtree_files, value = TRUE)
    tree_mixture_t_iqfile <- grep(paste0(a_m_prefix,".T"), all_iqtree_files, value = TRUE)
    
    # Extract information about each tree mixture model run
    tr_results <- extract.tree.mixture.results(tree_mixture_file = tree_mixture_tr_iqfile, 
                                               dataset = a_dataset, prefix = paste0(a_m_prefix,".TR"), 
                                               model = m, best_model = best_model, tree_branch_option = "TR")
    t_results <- extract.tree.mixture.results(tree_mixture_file = tree_mixture_t_iqfile, 
                                              dataset = a_dataset, prefix = paste0(a_m_prefix,".T"), 
                                              model = m,  best_model = best_model, tree_branch_option = "T")
    
    # Collate results dataframes
    op_df <- rbind(tr_results, r_results)
    names(op_df) <- names(tr_results)
    
    # Output results dataframe
    op_file <- paste0(a_tm_op_dir, a_m_prefix, "_tree_mixture_results.csv")
    write.csv(op_df, file = op_file, row.names = FALSE)
    
  } # end for (m in model_components)
} # end for (a in all_alignments)
