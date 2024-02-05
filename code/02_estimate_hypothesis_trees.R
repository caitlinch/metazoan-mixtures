## caitlinch/metazoan-mixtures/code/02_estimate_hypothesis_trees.R
# This script estimates maximum likelihood under constraint trees for 14 empirical data sets
# Caitlin Cherryh 2023


#### 1. Input parameters ####
## Specify parameters:
# alignment_dir               <- Directory containing alignments for all data sets
#                                   Alignment naming convention: [manuscript].[matrix_name].[sequence_type].fa
#                                   E.g. Cherryh2022.alignment1.aa.fa
# ml_tree_op_dir              <- Directory containing all maximum likelihood trees/log files/iqtree files estimated in pipeline step 01
# pmsf_op_dir                 <-  Directory containing all output files from estimating PMSF trees in pipeline step 01
# output_dir                  <- Directory for IQ-Tree output (trees and tree mixtures)
# repo_dir                    <- Location of caitlinch/metazoan-mixtures github repository
# number_parallel_processes   <- The number of simultaneous processes to run at once using mclapply(). 
#                                     If 1, then all processes will run sequentially
# iqtree2                     <- Location of IQ-Tree2 stable release
# iqtree_num_threads          <- Number of parallel threads for IQ-Tree to use. Can be a set number (e.g. 2) or "AUTO"
#                                     See IQ-Tree documentation for more details (http://www.iqtree.org/doc/Command-Reference)
# hypothesis_tree_bootstraps  <- Number of ultrafast bootstraps (UFB) to perform in IQ-Tree when estimating constrained maximum likelihood trees

## Specify control parameter values (all take logical values TRUE or FALSE):
# extract.ML.tree.information   <- TRUE to extract information from maximum likelihood tree log file and iqtree file, including tree topology. FALSE to skip.
# prepare.hypothesis.trees      <- TRUE to prepare constraint trees and create command lines to estimate hypothesis trees (constrained maximum likelihood trees). FALSE to skip.
# estimate.hypothesis.trees     <- TRUE to estimate all hypothesis trees (constrained maximum likelihood trees). FALSE to skip.

location = "local"
if (location == "local"){
  alignment_dir     <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/01_Data_all/"
  ml_tree_op_dir    <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/02_maximum_likelihood_trees/01_ml_tree_output_files/"
  pmsf_op_dir       <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/02_pmsf_site_freqs/"
  output_dir        <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/01_output_files/"
  repo_dir          <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"
  iqtree2           <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_Software_IQ-Tree/iqtree-2.2.0-MacOSX/bin/iqtree2"
  number_parallel_processes <- 1
  iqtree_num_threads        <- 3
} else if (location == "dayhoff"){
  alignment_dir     <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/data_all/"
  ml_tree_op_dir    <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/output/ml_tree_output_files/"
  pmsf_op_dir       <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/output/pmsf_trees/"
  output_dir        <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/output/"
  repo_dir          <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/"
  iqtree2           <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2"
  number_parallel_processes <- 4
  iqtree_num_threads        <- 20
}  else if (location == "rosa"){
  alignment_dir     <- "/home/caitlin/metazoan_mixtures/data_all/"
  ml_tree_op_dir    <- "/home/caitlin/metazoan_mixtures/output/ml_tree_output_files/"
  pmsf_op_dir       <- "/home/caitlin/metazoan_mixtures/output/pmsf_trees/"
  output_dir        <- "/home/caitlin/metazoan_mixtures/output/"
  repo_dir          <- "/home/caitlin/metazoan_mixtures/metazoan-mixtures/"
  iqtree2           <- "/home/caitlin/metazoan_mixtures/iqtree2/iqtree-2.2.2-Linux/bin/iqtree2"
  number_parallel_processes <- 4
  iqtree_num_threads        <- 20
}

# Create file paths for the constraint trees and hypothesis trees
c_tree_dir <- paste0(output_dir, "constraint_trees/")
h_tree_dir <- paste0(output_dir, "hypothesis_trees/")

# Set parameters that are identical for all run locations
hypothesis_tree_bootstraps <- NA

pmsf_model <- c("'LG+C60+F+R4'", "'LG+C20+F+R4'", "'C60+F+R4'", "'C20+F+R4'")
names(pmsf_model) <- c("PMSF_LG_C60", "PMSF_LG_C20", "PMSF_C60", "PMSF_C20")

# Set control parameters
control_parameters <- list(extract.ML.tree.information = FALSE,
                           prepare.constraint.trees = FALSE,
                           prepare.hypothesis.trees = FALSE,
                           estimate.hypothesis.trees = FALSE,
                           include.datasets.with.only.CXX.missing = TRUE)



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
   philippe2011_list, pick2010_list, ryan2013_list, simion2017_list, whelan2015_list, whelan2017_list, models_list, all_taxa, all_models)



#### 2. Prepare variables ####
# Extend the number of digits allowed (so BIC and logL can be extracted properly from iqtree files)
options(digits = 12)

# Create a folder for the constraint trees
if (file.exists(c_tree_dir) == FALSE){dir.create(c_tree_dir)}
# Create a folder for the hypothesis trees
if (file.exists(h_tree_dir) == FALSE){dir.create(h_tree_dir)}

# Create file paths for output files
output_file_paths <- paste0(output_dir, c("01_01_maximum_likelihood_tree_estimation_parameters.tsv",
                                          "01_01_all_tree_estimation_parameters.tsv",
                                          "01_02_maximum_likelihood_results.tsv",
                                          "01_02_maximum_likelihood_included_taxa.tsv",
                                          "01_02_dataset_completion_frequency.tsv",
                                          "01_03_best_models_per_alignment.tsv",
                                          "01_03_check_ModelFinder_best_models.tsv",
                                          "01_04_collated_constraint_tree_estimation_parameters.tsv",
                                          "01_04_constraint_tree_iqtree2_command_paths.txt"))



#### 4. Extract information from ML trees and log files ####
# Extract information about each run from the IQ-Tree output and log files
if (control_parameters$extract.ML.tree.information == TRUE){
  # Open ml_tree_df file 
  ml_tree_df <- read.table(output_file_paths[1], header = T)
  
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
  write.table(tree_df, file = output_file_paths[2], row.names = FALSE, sep = "\t")
  
  # Make a list of .iqtree files (and .log files)
  all_iqtree_files <- paste0(ml_tree_op_dir, tree_df$iqtree_file)
  all_log_files <- paste0(ml_tree_op_dir, gsub(".iqtree", ".log", tree_df$iqtree_file))
  
  # Determine which files exist (i.e. that ML tree estimated in IQ-Tree)
  finished_iqtree_files <- all_iqtree_files[file.exists(all_iqtree_files)]
  
  # Reduce the df to only rows that both the completed iqtree file and the log file are present for
  trimmed_ml_tree_df <- tree_df[tree_df$iqtree_file %in% basename(finished_iqtree_files),]
  rownames(trimmed_ml_tree_df) <- 1:nrow(trimmed_ml_tree_df)
  
  # Get the correct order for the .iqtree files by reading off the trimmed_ml_tree_df$iqtree_file column
  complete_iqtree_files <- paste0(ml_tree_op_dir, trimmed_ml_tree_df$iqtree_file)
  
  # Determine the completed log files - every finished tree will have a .treefile, a .iqtree file and a .log file
  complete_log_files <- paste0(ml_tree_op_dir, gsub(".iqtree", ".log", trimmed_ml_tree_df$iqtree_file))
  
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
  # Update best model for PMSF models
  pmsf_model <- c("'LG+C60+F+R4'", "'LG+C20+F+R4'", "'C60+F+R4'", "'C20+F+R4'")
  names(pmsf_model) <- c("PMSF_LG_C60", "PMSF_LG_C20", "PMSF_C60", "PMSF_C20")
  pmsf_rows <- which(grepl("PMSF", trimmed_ml_tree_df$model_code))
  trimmed_ml_tree_df$best_model[pmsf_rows] <- paste0(unlist(lapply(trimmed_ml_tree_df$model_code[pmsf_rows], function(x){pmsf_model[x]})),
                                                     ":", paste0(trimmed_ml_tree_df$prefix[pmsf_rows], ".ssfp.sitefreq"))
  
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
  # Extract the CXX parameters and format them nicely
  trimmed_ml_tree_df$estimated_CXX_frequencies <- unlist(lapply(complete_iqtree_files, extract.cat.frequencies, allow.zero.weights = TRUE))
  trimmed_ml_tree_df$estimated_CXX_frequencies_noZeroWeights <- unlist(lapply(complete_iqtree_files, extract.cat.frequencies, allow.zero.weights = FALSE))
  trimmed_ml_tree_df$estimated_CXX_frequencies_identical <- trimmed_ml_tree_df$estimated_CXX_frequencies == trimmed_ml_tree_df$estimated_CXX_frequencies_noZeroWeights
  
  # Update data frame to include maximum likelihood trees
  trimmed_ml_tree_df$maximum_likelihood_tree <- unlist(lapply(paste0(ml_tree_op_dir, trimmed_ml_tree_df$ml_tree_file), extract.treefile))
  
  # Reorder by dataset, then matrix name, then best tree by BIC
  trimmed_ml_tree_df <- trimmed_ml_tree_df[order(trimmed_ml_tree_df$dataset, trimmed_ml_tree_df$matrix_name, trimmed_ml_tree_df$tree_BIC),]
  
  # Add model class (Single, Other, PMSF, CXX)
  trimmed_ml_tree_df$model_class <- factor(trimmed_ml_tree_df$model_code,
                                           levels = c("PMSF_C60", "PMSF_LG_C60", "PMSF_C20", "PMSF_LG_C20", 
                                                      "LG_C60", "LG_C20", "C60", "C20", 
                                                      "CF4", "EHO", "EX2", "EX3", "EX_EHO", "LG4M", "UL2", "UL3", "ModelFinder",
                                                      "GTR20",  "JTTDCMut", "JTT", "LG", "mtZOA", "PMB", "Poisson", "rtREV", "WAG"),
                                           labels = c(rep("PMSF", 4), rep("CXX", 4), rep("Other", 9), rep("Single", 9) ))
  
  # Remove unwanted columns from the trimmed_ml_tree_df
  trimmed_ml_tree_df <- trimmed_ml_tree_df[,c("dataset", "model_code", "model_class", "matrix_name", "sequence_format",
                                              "prefix", "best_model", "best_model_LogL", "best_model_BIC", "best_model_wBIC",
                                              "tree_LogL", "tree_UnconstrainedLogL", "tree_NumFreeParams", "tree_BIC",
                                              "tree_length", "tree_SumInternalBranch", "tree_PercentInternalBranch",
                                              "estimated_rates", "estimated_gamma", "estimated_state_frequencies", 
                                              "estimated_CXX_frequencies", "estimated_CXX_frequencies_noZeroWeights",
                                              "estimated_CXX_frequencies_identical", "maximum_likelihood_tree")]
  # Save dataframe
  write.table(trimmed_ml_tree_df, file = output_file_paths[3], row.names = FALSE, sep = "\t")
  
  # Create a dataframe listing the taxa included in trees for each dataset
  if (file.exists(output_file_paths[4]) == FALSE){
    # Determine which taxa are included in the ML trees for each alignment (each value dataset/matrix name combination)
    alignment_taxa_df <- dataset.check.tree.taxa.wrapper(unique_ids = unique(paste0(trimmed_ml_tree_df$dataset, ".", trimmed_ml_tree_df$matrix_name)),
                                                         tree_folder = ml_tree_op_dir)
    # Save dataframe
    write.table(alignment_taxa_df, file = output_file_paths[4], row.names = FALSE, sep = "\t")
  }
}



#### 5. Estimate constraint and hypothesis trees for each combination of model and dataset ####
# Move to the folder for the constraint trees
setwd(c_tree_dir)

## Prepare constraint trees
if (control_parameters$prepare.constraint.trees == TRUE){
  ## Retrieve results from previous steps
  # Open output_file_paths[2] (input parameters to estimate trees in IQ-Tree)
  # Open trimmed_ml_tree_df file (output from ML tree runs)
  tree_df <- read.table(output_file_paths[2], header = T)
  # Open trimmed_ml_tree_df file (output from ML tree runs)
  trimmed_ml_tree_df <- read.table(output_file_paths[3], header = T)
  # Open alignment_taxa_df file (list of taxa in ML trees for each alignment)
  alignment_taxa_df <- read.table(output_file_paths[4], header = T)
  
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
  # Add another column
  completion_df$remaining_trees_to_run <- unlist(lapply(paste0(completion_df$dataset, ".", completion_df$matrix_name), check.remaining.runs, 
                                                        input_parameter_file = output_file_paths[2], output_parameter_file = output_file_paths[3]))[c(TRUE,FALSE)]
  completion_df$only_CXX_runs_remaining <- unlist(lapply(paste0(completion_df$dataset, ".", completion_df$matrix_name), check.remaining.runs, 
                                                         input_parameter_file = output_file_paths[2], output_parameter_file = output_file_paths[3]))[c(FALSE,TRUE)]
  # Output the frequency dataframe
  write.table(completion_df, file = output_file_paths[5], row.names = FALSE, sep = "\t")
  # Extract the names of the completed datasets/alignment combinations
  if (control_parameters$include.datasets.with.only.CXX.missing == TRUE){
    # Extract the names of the datasets/alignment combinations with EITHER:
    #       - all 26 models completed
    #       - all models EXCEPT any/all of the following completed: C20, C60, LG+C20, LG+C60
    complete_inds <- sort(unique(c(which(completion_df$frequency == 26), which(completion_df$only_CXX_runs_remaining == control_parameters$include.datasets.with.only.CXX.missing))))
  } else {
    # Extract the names of the datasets/alignment combinations with all 26 models completed
    complete_inds <- which(completion_df$frequency == 26)
  }
  completed_runs_df <- completion_df[complete_inds,]
  
  ## Determine which models to use for each completed dataset
  # Want to extract ModelFinder model, and the model with the best BIC
  # If the ModelFinder model has the best BIC, return just the ModelFinder model
  selected_models_df <- do.call(rbind, lapply(1:nrow(completed_runs_df), 
                                              determine.best.ML.model.per.class.wrapper, 
                                              completed_runs_df = completed_runs_df, 
                                              ML_output_df = trimmed_ml_tree_df))

  
  # Add the two substitution models to the csv
  # Nosenko 2013 non-ribo C60 has errors - replace with Nosenko 2013 non-ribo LG+C60
  selected_models_df <- rbind(selected_models_df, trimmed_ml_tree_df[which(trimmed_ml_tree_df$dataset == "Nosenko2013" & 
                                                                             trimmed_ml_tree_df$matrix_name == "nonribosomal_9187_smatrix" & 
                                                                             trimmed_ml_tree_df$model_code == "LG_C60"), ])
  selected_models_df <- rbind(selected_models_df, trimmed_ml_tree_df[which(trimmed_ml_tree_df$dataset == "Laumer2018" & 
                                                                             trimmed_ml_tree_df$matrix_name == "Tplx_phylo_d1" & 
                                                                             trimmed_ml_tree_df$model_code == "LG_C20"), ])
  
  # Save the dataframe of best models
  write.table(selected_models_df, file = output_file_paths[6], row.names = FALSE, sep = "\t")
  
  ## Construct constraint trees
  # Constraint and hypothesis trees will only be estimated for the best model(s) for each dataset/matrix combination (found in the selected_models_df)
  # Create the constraint trees and determine what parameters to use for each hypothesis tree
  #     Hypothesis tree = constrained maximum likelihood tree estimated in IQ-Tree with best model from ML run (for each constraint tree for each dataset/model combination)
  # Create a constraint dataframe for each row in the selected_models_df
  lapply(1:nrow(selected_models_df), output.all.constraint.trees, selected_models_df, 
         constraint_tree_dir = c_tree_dir, dataset_info = all_datasets, 
         matrix_taxa_info = matrix_taxa, ml_tree_tips_df = alignment_taxa_df)
}

## Prepare hypothesis tree commands
if (prepare.hypothesis.trees == TRUE){
  ### Prepare hypothesis tree estimation
  # Copy df for hypothesis tree estimation
  hyp_te_df <- read.table(output_file_paths[6], header = TRUE, sep = "\t")
  # Attach alignment files to the rows
  all_alignment_files <- list.files(alignment_dir)
  all_alignment_files <- grep("alignment", grep("00_", all_alignment_files, invert = T, value = T), value = T)
  hyp_te_df$alignment_file <- unlist(lapply(1:nrow(hyp_te_df), function(i){
    grep(hyp_te_df$dataset[i], grep(hyp_te_df$matrix_name[i], all_alignment_files, value = T), value = T)}))
  # Add the mrate = NA options for IQ-Tree to the dataframe (do not include mrate option for estimating constraint trees)
  hyp_te_df$model_mrate <- NA
  # Add the number of threads
  hyp_te_df$num_threads <- iqtree_num_threads
  # Add the number of ultrafast bootstraps to perform for the hypothesis trees (constrained maximum likelihood trees with fixed model)
  hyp_te_df$num_bootstraps <- hypothesis_tree_bootstraps
  # Add new column with the number of constraint trees to estimate
  hyp_te_df$num_constraint_trees <- 5
  hyp_te_df[hyp_te_df$dataset == "Borowiec2015", ]$num_constraint_trees <- 3
  hyp_te_df[hyp_te_df$dataset == "Dunn2008", ]$num_constraint_trees <- 3
  # Add column with constraint path name by duplicating each row the number of times required (once per tree)
  c_te_df <- duplicate.constraint.rows(df = hyp_te_df)
  c_te_df$constraint_tree_path <- paste0(c_te_df$dataset, ".", c_te_df$matrix_name, ".constraint_tree_", c_te_df$constraint_tree_id, ".nex")
  c_te_df$hypothesis_tree_path <- paste0(c_te_df$dataset, ".", c_te_df$matrix_name, ".",  c_te_df$model_code, ".ML_H", c_te_df$constraint_tree_id)
  # Update filepaths for the constraint_df
  c_te_df$alignment_file        <- paste0(alignment_dir, basename(c_te_df$alignment_file))
  c_te_df$constraint_tree_path <- paste0(c_tree_dir, basename(c_te_df$constraint_tree_path))
  c_te_df$hypothesis_tree_path <- paste0(c_tree_dir, basename(c_te_df$hypothesis_tree_path))
  c_te_df$iqtree_path <- iqtree2
  # Update best_model column for PMSF models to add file path to site frequency files
  pmsf_model <- c("'LG+C60+F+R4'", "'LG+C20+F+R4'", "'C60+F+R4'", "'C20+F+R4'")
  names(pmsf_model) <- c("PMSF_LG_C60", "PMSF_LG_C20", "PMSF_C60", "PMSF_C20")
  pmsf_rows <- which(grepl("PMSF", c_te_df$model_code))
  c_te_df$best_model[pmsf_rows] <- paste0(unlist(lapply(c_te_df$model_code[pmsf_rows], function(x){pmsf_model[x]})), 
                                          ":", paste0(pmsf_op_dir, c_te_df$prefix[pmsf_rows], ".ssfp.sitefreq"))
  
  # Prepare iqtree commands for each of the hypothesis trees
  c_te_df$iqtree2_call <- unlist(lapply(1:nrow(c_te_df), construct.hypothesis.tree.call, hyp_tree_info_df = c_te_df))
  
  # Save the constraint tree dataframe
  write.table(c_te_df, file = output_file_paths[8], row.names = FALSE, sep = "\t")
  # Save list of iqtree2 commands as text file
  write(c_te_df$iqtree2_call, file = output_file_paths[9])
}

## Estimate hypothesis trees
if (control_parameters$estimate.hypothesis.trees == TRUE){
  # Open constraint tree dataframe file
  constraint_df <- read.table(output_file_paths[8], header = T)
  
  # Run IQ-Tree commands to estimate hypothesis trees for each model/matrix combination
  mclapply(constraint_df$iqtree2_call, system, mc.cores = number_parallel_processes)
}



#### 6. Collate hypothesis trees into one file per dataset/model combination ####
# List all files in the hypothesis tree directory
all_files <- list.files(h_tree_dir, recursive = T)
# Remove collated hypothesis trees
all_files <- grep("collated_hypothesis_trees|hypothesis|constrained_ML", all_files, value = T, invert = T)
# Get the files for hypothesis trees
filenames <- basename(all_files)
iq_files <- grep("\\.iqtree", filenames, value = T) # get all IQ-Tree files
ml_h_files <- grep("ML_H1", iq_files, value = T) # get only files for hypothesis trees
# Separate out the unique IDs for each combination of model and dataset
run_ids <- sort(unique(unlist(lapply(strsplit(ml_h_files, "\\."), function(x){paste0(x[[1]], ".", x[[2]], ".", x[[3]])}))))
# Collate hypothesis trees by id
collated_h_2trees <- lapply(run_ids, combine.hypothesis.trees, tree_directory = h_tree_dir, output_id = "hypothesis_trees.2trees", num_trees = 2)
collated_h_5trees <- lapply(run_ids, combine.hypothesis.trees, tree_directory = h_tree_dir, output_id = "hypothesis_trees.5trees", num_trees = 5)


