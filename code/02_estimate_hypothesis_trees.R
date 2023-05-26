# metazoan-mixtures/code/02_estimate_hypothesis_trees.R
## This script estimates maximum likelihood under constraint trees for 14 empirical data sets
# Caitlin Cherryh, 2022



#### 1. Input parameters ####
## Specify parameters:
# alignment_dir               <- Directory containing alignments for all data sets
#                                   Alignment naming convention: [manuscript].[matrix_name].[sequence_type].fa
#                                   E.g. Cherryh2022.alignment1.aa.fa
# ml_tree_dir                 <- Directory containing all maximum likelihood trees estimating in pipeline step 01
# output_dir                  <- Directory for IQ-Tree output (trees and tree mixtures)
# repo_dir                    <- Location of caitlinch/metazoan-mixtures github repository
# number_parallel_processes   <- The number of simultaneous processes to run at once using mclapply(). 
#                                     If 1, then all processes will run sequentially
# iqtree2                     <- Location of IQ-Tree2 stable release
# iqtree_num_threads          <- Number of parallel threads for IQ-Tree to use. Can be a set number (e.g. 2) or "AUTO"
#                                     See IQ-Tree documentation for more details (http://www.iqtree.org/doc/Command-Reference)
# hypothesis_tree_bootstraps  <- Number of ultrafast bootstraps (UFB) to perform in IQ-Tree when estimating constrained maximum likelihood trees

## Specify control parameters (all take logical values TRUE or FALSE:
# extract.ML.tree.information   <- TRUE to extract information from maximum likelihood tree log file and iqtree file, including tree topology. FALSE to skip.
# prepare.hypothesis.trees      <- TRUE to prepare constraint trees and create command lines to estimate hypothesis trees (constrained maximum likelihood trees). FALSE to skip.
# estimate.hypothesis.trees     <- TRUE to estimate all hypothesis trees (constrained maximum likelihood trees). FALSE to skip.
# collate.hypothesis.logs       <- TRUE to extract information from hypothesis tree log file and iqtree file. FALSE to skip.

location = "local"
if (location == "local"){
  alignment_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/01_Data_all/"
  ml_tree_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/02_maximum_likelihood_trees/01_ml_tree_output_files/"
  output_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/"
  repo_dir <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"
  number_parallel_processes <- 1
  iqtree2 <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_Software_IQ-Tree/iqtree-2.2.0-MacOSX/bin/iqtree2"
} else if (location == "dayhoff"){
  alignment_dir <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/data_all/"
  ml_tree_dir <- ""
  output_dir <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/output/"
  repo_dir <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/"
  number_parallel_processes <- 4
  iqtree2 <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2"
} 

# Set parameters that are identical for all run locations
iqtree_num_threads <- 15
hypothesis_tree_bootstraps <- 1000

# Set control parameters
extract.ML.tree.information <- FALSE
prepare.hypothesis.trees <- FALSE
estimate.hypothesis.trees <- FALSE
collate.hypothesis.logs <- FALSE



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
# Create output folders
# Create a folder for the constraint trees
c_tree_dir <- paste0(output_dir, "constraint_trees/")
if (file.exists(c_tree_dir) == FALSE){dir.create(c_tree_dir)}
# Create a folder for the hypothesis trees
h_tree_dir <- paste0(output_dir, "hypothesis_trees/")
if (file.exists(h_tree_dir) == FALSE){dir.create(h_tree_dir)}

# Create file paths for output files
ml_tree_df_file             <- paste0(repo_dir, "output/01_01_maximum_likelihood_tree_estimation_parameters.tsv")
all_tree_df_file            <- paste0(repo_dir, "output/01_01_all_tree_estimation_parameters.tsv")
ml_extracted_df_file        <- paste0(repo_dir, "output/01_02_maximum_likelihood_results.tsv")
alignment_taxa_df_file      <- paste0(repo_dir, "01_02_maximum_likelihood_included_taxa.tsv")
completion_freq_df_file     <- paste0(repo_dir, "01_02_dataset_completion_frequency.tsv")



#### 4. Extract information from ML trees and log files ####
# Extract information about each run from the IQ-Tree output and log files
if (extract.ML.tree.information == TRUE){
  # Open ml_tree_df file 
  ml_tree_df <- read.table(ml_tree_df_file, header = T)
  
  # Add the PMSF runs
  pmsf_ids <- grep("LG_C20|LG_C60|C20|C60", ml_tree_df$model_code)
  pmsf_df <- ml_tree_df[pmsf_ids,]
  pmsf_df$iqtree2_call <- NA
  pmsf_df$model_mrate <- NA
  pmsf_df$model_code <- paste0("PMSF_", pmsf_df$model_code)
  pmsf_df$model_mset <- paste0("PMSF ", pmsf_df$model_mset)
  pmsf_df$prefix <- paste0(pmsf_df$dataset, ".", pmsf_df$matrix_name, ".", pmsf_df$model_code)
  pmsf_df$iqtree_file <- paste0(pmsf_df$prefix, ".complete.iqtree")
  pmsf_df$ml_tree_file <- paste0(pmsf_df$prefix, ".complete.treefile")
  
  # Attach pmsf and ml tree runs into one big dataset
  tree_df <- rbind(ml_tree_df, pmsf_df)
  # Sort by dataset then by model
  tree_df <- tree_df[order(tree_df$dataset, tree_df$matrix_name, tree_df$model_code),]
  # Write the combined table out
  write.table(tree_df, file = all_tree_df_file, row.names = FALSE, sep = "\t")
  
  # Make a list of .iqtree files (and .log files)
  all_iqtree_files <- paste0(ml_tree_dir, tree_df$iqtree_file)
  all_log_files <- paste0(ml_tree_dir, gsub(".iqtree", ".log", tree_df$iqtree_file))
  
  # Determine which files exist (i.e. that ML tree estimated in IQ-Tree)
  finished_iqtree_files <- all_iqtree_files[file.exists(all_iqtree_files)]
  
  # Reduce the df to only rows that both the completed iqtree file and the log file are present for
  trimmed_ml_tree_df <- tree_df[tree_df$iqtree_file %in% basename(finished_iqtree_files),]
  rownames(trimmed_ml_tree_df) <- 1:nrow(trimmed_ml_tree_df)
  
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
  trimmed_ml_tree_df$tree_PercentInternalBranch <- round(((as.numeric(trimmed_ml_tree_df$tree_SumInternalBranch)/as.numeric(trimmed_ml_tree_df$tree_length)) * 100), digits = 2)
  
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
  
  # Remove unwanted columns from the trimmed_ml_tree_df
  trimmed_ml_tree_df <- trimmed_ml_tree_df[,c("dataset", "model_code", "matrix_name", "sequence_format", "prefix", 
                                              "tree_LogL", "tree_UnconstrainedLogL", "tree_NumFreeParams", "tree_BIC",
                                              "tree_length", "tree_SumInternalBranch", "tree_PercentInternalBranch",
                                              "best_model", "best_model_LogL", "best_model_BIC", "best_model_wBIC",
                                              "estimated_rates", "estimated_gamma", "estimated_state_frequencies", 
                                              "maximum_likelihood_tree")]
  # Save dataframe
  write.table(trimmed_ml_tree_df, file = ml_extracted_df_file, row.names = FALSE, sep = "\t")
  
  # Create a dataframe listing the taxa included in trees for each dataset
  if (file.exists(alignment_taxa_df_file) == FALSE){
    # Determine which taxa are included in the ML trees for each alignment (each value dataset/matrix name combination)
    alignment_taxa_df <- dataset.check.tree.taxa.wrapper(unique_ids = unique(paste0(trimmed_ml_tree_df$dataset, ".", trimmed_ml_tree_df$matrix_name)),
                                                         tree_folder = ml_tree_dir)
    # Save dataframe
    write.table(alignment_taxa_df, file = alignment_taxa_df_file, row.names = FALSE, sep = "\t")
  }
}



#### 5. Estimate constraint and hypothesis trees for each combination of model and dataset ####
# Move to the folder for the constraint trees
setwd(c_tree_dir)

# Prepare constraint trees to estimate hypothesis trees
if (prepare.hypothesis.trees == TRUE){
  ## Retrieve results from previous steps
  # Open trimmed_ml_tree_df file (output from ML tree runs)
  trimmed_ml_tree_df <- read.table(ml_extracted_df_file, header = T)
  # Open alignment_taxa_df file (list of taxa in ML trees for each alignment)
  alignment_taxa_df <- read.table(alignment_taxa_df_file, header = T)
  
  ## Select completed datasets to estimate constraint trees
  # Determine which datasets have all alignments completed
  completion_df <- as.data.frame(table(trimmed_ml_tree_df$dataset, trimmed_ml_tree_df$matrix_name), stringsAsFactors = FALSE)
  names(completion_df) <- c("dataset", "matrix_name", "frequency")
  # Remove all entries with 0 frequency (either an alignment that was not run, or an artifact of the method for making the table 
  #   i.e. a combination of dataset and alignment name that is incorrect)
  completion_df <- completion_df[completion_df$frequency != 0,]
  # Sort by dataset then by matrix
  completion_df <- completion_df[order(completion_df$dataset, completion_df$matrix_name),]
  # Reset row names 
  row.names(completion_df) <- 1:nrow(completion_df)
  # Output the frequency dataframe
  write.table(completion_df, file = completion_freq_df_file, row.names = FALSE, sep = "\t")
  # Extract the names of the datasets/alignment combinations with all 24 models completed
  completed_df <- completion_df[completion_df$frequency == 26,]
  
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


