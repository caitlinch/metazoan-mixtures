# metazoan-mixtures/code/02_estimate_hypothesis_trees.R
## This script estimates maximum likelihood under constraint trees for 14 empirical data sets
# Caitlin Cherryh, 2022

## This script estimates ML trees using a constraint tree for empirical data sets under different models



#### 1. Input parameters ####
## Specify parameters:
# alignment_dir       <- Directory containing alignments for all data sets
#                           Alignment naming convention: [manuscript].[matrix_name].[sequence_type].fa
#                           E.g. Cherryh2022.alignment1.aa.fa
# output_dir          <- Directory for IQ-Tree output (trees and tree mixtures)
# repo_dir            <- Location of caitlinch/metazoan-mixtures github repository

# iqtree2             <- Location of IQ-Tree2 stable release

# iqtree_num_threads          <- Number of parallel threads for IQ-Tree to use. Can be a set number (e.g. 2) or "AUTO"
# iqtree_mrate                <- Specify a comma separated list of rate heterogeneity types for model selection in IQ-Tree
#                                   We set iqtree_mrate = "E,I,G,I+G,R,I+R"
#                                   See IQ-Tree documentation for more details (http://www.iqtree.org/doc/Command-Reference)
# ml_tree_bootstraps          <- Number of ultrafast bootstraps (UFB) to perform in IQ-Tree
# hypothesis_tree_bootstraps  <- Number of ultrafast bootstraps (UFB) to perform in IQ-Tree when estimating constrained maximum likelihood trees
# number_parallel_processes   <- The number of simultaneous processes to run at once using mclapply(). 
#                                   If 1, then all processes will run sequentially


## Specify control parameters (all take logical values TRUE or FALSE:
# estimate.ML.trees             <- TRUE to estimate all maximum likelihood trees (each combination of model and alignment). FALSE to skip.
# extract.ML.tree.information   <- TRUE to extract information from maximum likelihood tree log file and iqtree file, including tree topology. FALSE to skip.
# prepare.hypothesis.trees      <- TRUE to prepare constraint trees and create command lines to estimate hypothesis trees (constrained maximum likelihood trees). FALSE to skip.
# estimate.hypothesis.trees     <- TRUE to estimate all hypothesis trees (constrained maximum likelihood trees). FALSE to skip.
# collate.hypothesis.logs       <- TRUE to extract information from hypothesis tree log file and iqtree file. FALSE to skip.

location = "local"
if (location == "local"){
  ml_tree_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/01_Data_all/"
  output_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/"
  repo_dir <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"
  
  iqtree2 <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_Software_IQ-Tree/iqtree-2.2.0-MacOSX/bin/iqtree2"
  
  number_parallel_processes <- 1
  
} else if (location == "soma"){
  ml_tree_dir <- ""
  output_dir <- "/data/caitlin/metazoan-mixtures/output/"
  repo_dir <- "/data/caitlin/metazoan-mixtures/"
  
  iqtree2 <- "/data/caitlin/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2"
  
  number_parallel_processes <- 2
  
} else if (location == "dayhoff"){
  ml_tree_dir <- ""
  output_dir <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/output/"
  repo_dir <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/"
  
  iqtree2 <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2"
  
  number_parallel_processes <- 1
  
} else if (location == "laptop"){
  ml_tree_dir <- ""
  output_dir <- "/Users/caitlin/Documents/PhD/Ch03_sponge_mixtures/02_output/"
  repo_dir <- "/Users/caitlin/Repositories/metazoan-mixtures/"
  
  iqtree2 <- "/Users/caitlin/Documents/PhD/Ch03_sponge_mixtures/iqtree-2.2.0-MacOSX/bin/iqtree2"
  
  number_parallel_processes <- 1
} else if (location == "rona"){
  ml_tree_dir <- ""
  output_dir <- "/home/caitlin/metazoan-mixtures/output/"
  repo_dir <- "/home/caitlin/metazoan-mixtures/"
  
  iqtree2 <- "/home/caitlin/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2"
  
  number_parallel_processes <- 1
}

# Set parameters that are identical for all run locations
iqtree_mrate <- "E,I,G,I+G,R,I+R"
iqtree_num_threads <- 15
ml_tree_bootstraps <- 1000
hypothesis_tree_bootstraps <- 1000

# Set control parameters
estimate.ML.trees <- FALSE
extract.ML.tree.information <- FALSE
prepare.hypothesis.trees <- FALSE # Completed this step for Borowiec2015, Nosenko2013 nonribo, Philippe2009
estimate.hypothesis.trees <- TRUE
collate.hypothesis.logs <- TRUE



#### 2. Prepare functions and packages ####
# Open packages
library(parallel)
library(ape)

# Source functions
source(paste0(repo_dir, "code/func_estimate_trees.R"))
source(paste0(repo_dir, "code/func_data_processing.R"))

# Source information about datasets
source(paste0(repo_dir, "code/data_dataset_info.R"))

# Remove the individual dataset lists (only need collated lists) (yes it is a bit cheeky to hard code the removal)
rm(borowiec2015_list, chang2015_list, dunn2008_list, hejnol2009_list, laumer2018_list, laumer2019_list, moroz2014_list, nosenko2013_list, philippe2009_list,
   philippe2011_list, pick2010_list, ryan2013_list, simion2017_list, whelan2015_list, whelan2017_list, models_list, all_taxa)



#### 2. Prepare variables ####
# Extract the list of all files from the folder containing alignments/models/partition files
all_files <- list.files(alignment_dir)
if (length(all_files) > 0){
  all_files <- paste0(alignment_dir, all_files)
}
# Extract the list of alignments (i.e. files that contain the word "alignment")
all_alignments <- grep("\\.alignment\\.", all_files, value = T)

# Sort and categorise models of sequence evolution 
all_models <- sort(unique(unlist(lapply(all_models, sort.model.chunks)))) 
# Identify unique matrices within the models
model_components <- sort(unique(unlist(strsplit(all_models, "\\+"))))
# Remove components that do not correspond to a matrix or model (e.g. remove "I", "F", "G4")
model_components <- model_components[!model_components %in% c("F", "FO", "G", "G4", "I", "R", "R4")]
# Create a vector of models to exclude from analysis
#   PMSF: PMSF is more complex to run than other models - requires separate function/multiple steps (see separate PMSF file)
#   CAT: CAT is not a model in IQ-Tree. Use C10:C60 instead.
#   GTR: GTR is a DNA model
#   F81: F81 is a DNA model
exclude_models <- c("PMSF", "CAT", "GTR", "F81", "C10", "C30", "C40", "C50")
# Remove the models to exclude from the list of models
model_components <- model_components[which(!model_components %in% exclude_models)]
# Add the newer models 
model_components <- c(model_components, "'LG+C20'", "'LG+C60'")
# Sort the mode components
model_components <- sort(model_components)
# Move ModelFinder to the end of the list (will take the longest as has to test all possible models, so run last)
model_components <- model_components[!model_components == "ModelFinder"]
model_components <- c(model_components, "ModelFinder")
# Note: partitioning schemes currently not possible in mixture of trees implementation

# Create output folders
# Create a folder for the ml trees
ml_tree_dir <- paste0(output_dir, "maximum_likelihood_trees/")
if (file.exists(ml_tree_dir) == FALSE){dir.create(ml_tree_dir)}
# Create a folder for the constraint trees
c_tree_dir <- paste0(output_dir, "constraint_trees/")
if (file.exists(c_tree_dir) == FALSE){dir.create(c_tree_dir)}

# Create file paths for output files
ml_run_output_path <- paste0(output_dir, "01_02_maximum_likelihood_results.tsv")



#### 4. Extract information from ML trees and log files ####
# Extract information about each run from the IQ-Tree output and log files
if (extract.ML.tree.information == TRUE){
  # Open ml_tree_df file 
  ml_tree_df <- read.table(df_op_01_01, header = T)
  
  # Make a list of .iqtree files (and .log files)
  all_iqtree_files <- paste0(ml_tree_dir, ml_tree_df$iqtree_file)
  all_log_files <- paste0(ml_tree_dir, gsub(".iqtree", ".log", ml_tree_df$iqtree_file))
  
  # Determine which files exist (i.e. that ML tree estimated in IQ-Tree)
  finished_iqtree_files <- all_iqtree_files[file.exists(all_iqtree_files)]
  
  # Reduce the df to only rows that both the completed iqtree file and the log file are present for
  trimmed_ml_tree_df <- ml_tree_df[ml_tree_df$iqtree_file %in% basename(finished_iqtree_files),]
  
  # Get the correct order for the .iqtree files by reading off the trimmed_ml_tree_df$iqtree_file column
  complete_iqtree_files <- paste0(ml_tree_dir, trimmed_ml_tree_df$iqtree_file)
  
  # Determine the completed log files - every finished tree will have a .treefile, a .iqtree file and a .log file
  complete_log_files <- paste0(ml_tree_dir, gsub(".iqtree", ".log", trimmed_ml_tree_df$iqtree_file))
  
  # Extract the log likelihood and other values for the tree
  trimmed_ml_tree_df$tree_LogL <- unlist(lapply(complete_iqtree_files, extract.tree.log.likelihood, var = "LogL"))
  trimmed_ml_tree_df$tree_UnconstrainedLogL <- unlist(lapply(complete_iqtree_files, extract.tree.log.likelihood, var = "ULL"))
  trimmed_ml_tree_df$tree_NumFreeParams <- unlist(lapply(complete_iqtree_files, extract.tree.log.likelihood, var = "NFP"))
  trimmed_ml_tree_df$tree_BIC <- unlist(lapply(complete_iqtree_files, extract.tree.log.likelihood, var = "BIC"))
  trimmed_ml_tree_df$tree_length <- unlist(lapply(complete_iqtree_files, extract.tree.log.likelihood, var = "TrLen"))
  trimmed_ml_tree_df$tree_SumInternalBranch <- unlist(lapply(complete_iqtree_files, extract.tree.log.likelihood, var = "SIBL"))
  
  # Extract the best model for each combination of matrix and model
  trimmed_ml_tree_df$best_model <- unlist(lapply(complete_iqtree_files, extract.best.model))
  # Update the best model column to remove "+I+I+" values (should be "+I+" - output error from IQ-Tree)
  trimmed_ml_tree_df$best_model <- gsub("\\+I\\+I\\+", "+I+", trimmed_ml_tree_df$best_model)
  
  # Extract the BIC value and log likelihood value for the best model
  trimmed_ml_tree_df$best_model_LogL <- unlist(lapply(complete_iqtree_files, extract.model.log.likelihood, var = "LogL"))
  trimmed_ml_tree_df$best_model_BIC <- unlist(lapply(complete_iqtree_files, extract.model.log.likelihood, var = "BIC"))
  trimmed_ml_tree_df$best_model_wBIC <- unlist(lapply(complete_iqtree_files, extract.model.log.likelihood, var = "wBIC"))
  
  # Extract details about the model
  trimmed_ml_tree_df$estimated_rates <- unlist(lapply(complete_iqtree_files, extract.rates))
  trimmed_ml_tree_df$estimated_gamma <- unlist(lapply(complete_iqtree_files, extract.gamma.values))
  trimmed_ml_tree_df$estimated_state_frequencies <- unlist(lapply(complete_iqtree_files, extract.state.frequencies))
  
  # Update data frame to include maximum likelihood trees
  trimmed_ml_tree_df$maximum_likelihood_tree <- unlist(lapply(paste0(ml_tree_dir, trimmed_ml_tree_df$ml_tree_file), extract.treefile))
  
  # Save dataframe
  write.table(trimmed_ml_tree_df, file = df_op_01_02, row.names = FALSE, sep = "\t")
  
  # Determine which taxa are included in the ML trees for each alignment (each value dataset/matrix name combination)
  alignment_taxa_df <- dataset.check.tree.taxa.wrapper(unique_ids = unique(paste0(trimmed_ml_tree_df$dataset, ".", trimmed_ml_tree_df$matrix_name)),
                                                       tree_folder = ml_tree_dir)
  # Save dataframe
  write.table(alignment_taxa_df, file = df_op_alignment_taxa, row.names = FALSE, sep = "\t")
}



#### 5. Estimate constraint and hypothesis trees for each combination of model and dataset ####
# Move to the folder for the constraint trees
setwd(c_tree_dir)

# Prepare constraint trees to estimate hypothesis trees
if (prepare.hypothesis.trees == TRUE){
  ## Retrieve results from previous steps
  # Open trimmed_ml_tree_df file (output from ML tree runs)
  trimmed_ml_tree_df <- read.table(df_op_01_02, header = T)
  # Open alignment_taxa_df file (list of taxa in ML trees for each alignment)
  alignment_taxa_df <- read.table(df_op_alignment_taxa, header = T)
  
  ## Select completed datasets to estimate constraint trees
  # Determine which datasets have all alignments completed
  completion_df <- as.data.frame(table(trimmed_ml_tree_df$dataset, trimmed_ml_tree_df$matrix_name), stringsAsFactors = FALSE)
  names(completion_df) <- c("dataset", "matrix_name", "frequency")
  # Remove all entries with 0 frequency (either an alignment that was not run, or an artifact of the method for making the table 
  #   i.e. a combination of dataset and alignment name that is incorrect)
  completion_df <- completion_df[completion_df$frequency != 0,]
  row.names(completion_df) <- 1:nrow(completion_df)
  # Output the frequency dataframe
  write.table(completion_df, file = df_op_completion_freq, row.names = FALSE, sep = "\t")
  # Extract the names of the datasets/alignment combinations with all 24 models completed
  completed_df <- completion_df[completion_df$frequency == 24,]
  
  ## Determine which models to use for each completed dataset
  # Want to extract ModelFinder model, and the model with the best BIC
  # If the ModelFinder model has the best BIC, return just the ModelFinder model
  selected_models_list <- lapply(1:nrow(completed_df), determine.best.ML.model.wrapper, completed_runs_df = completed_df, 
                                 ML_output_df = trimmed_ml_tree_df, include.ModelFinder = FALSE) 
  # Convert lapply output to a nice dataframe
  selected_models_df <- do.call(rbind, selected_models_list)
  # Save the dataframe of best models
  write.table(selected_models_df, file = df_op_best_models, row.names = FALSE, sep = "\t")
  
  ## Check whether the "best model" by BIC is tested for by ModelFinder in IQ-Tree
  check_modelfinder_df <- selected_models_list <- do.call(rbind, lapply(1:nrow(completed_df), determine.best.ML.model.wrapper, completed_runs_df = completed_df, 
                                                                        ML_output_df = trimmed_ml_tree_df, include.ModelFinder = TRUE))
  mfp_check_df <- check.ModelFinder.models.wrapper(best_models_df = check_modelfinder_df, IQTree_output_dir = ml_tree_dir)
  # Save the dataframe you just created
  write.table(mfp_check_df, file = df_op_mfp_model_check, row.names = FALSE, sep = "\t")
  
  ## Prepare parameters for the constraint trees
  # Constraint and hypothesis trees will only be estimated for the best model(s) for each dataset/matrix combination (found in the selected_models_df)
  # Create the constraint trees and determine what parameters to use for each hypothesis tree
  #     Hypothesis tree = constrained maximum likelihood tree estimated in IQ-Tree with best model from ML run (for each constraint tree for each dataset/model combination)
  # Create a constraint dataframe for each row in the selected_models_df
  constraint_list <- lapply(1:nrow(selected_models_df), constraint.tree.wrapper, output_directory = c_tree_dir,
                            iqtree_path = iqtree2, iqtree_num_threads = iqtree_num_threads,
                            dataset_info = all_datasets, matrix_taxa_info = matrix_taxa,
                            ml_output_df = selected_models_df, ml_tree_tips_df = alignment_taxa_df, 
                            force.update.constraint.trees = TRUE)
  # Collate the constraints into a single dataframe
  constraint_df <- do.call(rbind, constraint_list)
  # Add the mrate = NA options for IQ-Tree to the dataframe (do not include mrate option for estimating constraint trees)
  constraint_df$model_mrate <- NA
  # Add the number of ultrafast bootstraps to perform for the hypothesis trees (constrained maximum likelihood trees with fixed model)
  constraint_df$num_bootstraps <- hypothesis_tree_bootstraps
  
  # Prepare iqtree commands for each of the hypothesis trees
  constraint_df$iqtree2_call <- unlist(lapply(1:nrow(constraint_df), run.one.constraint.tree, constraint_df = constraint_df, run.iqtree = FALSE))
  
  # Save the constraint tree dataframe
  write.table(constraint_df, file = df_op_01_03, row.names = FALSE, sep = "\t")
  
  # Save list of iqtree2 commands as text file
  write(constraint_df$iqtree2_call, file = txt_op_01_03)
}

# Estimate hypothesis trees
if (estimate.hypothesis.trees == TRUE){
  # Open constraint tree dataframe file
  constraint_df <- read.table(df_op_01_03, header = T)
  
  # Run IQ-Tree commands to estimate hypothesis trees for each model/matrix combination
  mclapply(constraint_df$iqtree2_call, system, mc.cores = number_parallel_processes)
}

# Extract hypothesis trees from output and save into a tsv
#### THIS NEEDS TESTING ####
if (collate.hypothesis.logs == TRUE){
  # Open ml_tree_df file (if it exists)
  if (file.exists(df_op_01_01) == TRUE){
    ml_tree_df <- read.table(df_op_01_01, header = T)
  }
  
  # Combine hypothesis trees into one file per unique identifier and save
  ml_tree_df$hypothesis_tree_files <- lapply(ml_tree_df$prefix, combine.hypothesis.trees, constraint_tree_directory = c_tree_dir,
                                             outgroup_taxa = NA)
  # Save dataframe
  write.table(ml_tree_df, file = df_op_01_04, row.names = FALSE, sep = "\t")
  
}


