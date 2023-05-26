# metazoan-mixtures/code/01_estimate_PMSF_trees.R
## This script estimates maximum likelihood trees under the PMSF model for 14 empirical data sets
# Caitlin Cherryh, 2022

## This script:
# 1. Estimates ML trees for empirical data sets under the PMSF model (and eventually under the CAT-PMSF models)
# 2. Estimates ML trees using a constraint tree for empirical data sets under the PMSF model (and eventually under the CAT-PMSF models)

# When using the posterior mean site frequency model (PMSF model), cite the following:
#   H.C. Wang, B.Q. Minh, S. Susko and A.J. Roger (2018), 
#     Modeling site heterogeneity with posterior mean site frequency profiles
#     accelerates accurate phylogenomic estimation. Syst. Biol., 67:216-235.
#     https://doi.org/10.1093/sysbio/syx068



#### 1. Input parameters ####
## Specify parameters:
# alignment_dir               <- Directory containing alignments for all data sets
#                                   Alignment naming convention: [manuscript].[matrix_name].[sequence_type].fa
#                                   E.g. Cherryh2022.alignment1.aa.fa
# output_dir                  <- Directory for IQ-Tree output (trees and tree mixtures)
# repo_dir                    <- Location of caitlinch/metazoan-mixtures github repository
# iqtree2                     <- Location of IQ-Tree2 stable release

# pmsf_initial_model          <- Model used to estimate guide tree for the site-specific frequency model (we use pmsf_initial_model = "'LG+C60+F+G'")
#                                   Place model inside two different types of quotes so it can be pasted into an IQ-Tree command line properly, e.g. "'model+R6'"
# iqtree_num_threads          <- Number of parallel threads for IQ-Tree to use. Can be a set number (e.g. 2) or "AUTO"
# ml_tree_bootstraps          <- Number of ultrafast bootstraps (UFB) to perform in IQ-Tree
# number_parallel_processes   <- The number of simultaneous processes to run at once using mclapply(). 
#                                   If 1, then all processes will run sequentially

# run_IQTREE                  <- FALSE to output IQ-Tree command lines only, TRUE to output IQ-Tree command lines and run IQ-Tree 
# trim_incomplete             <- TRUE to remove dataframe rows without a site frequency file before estimating tree in IQ-Tree
#                                   (tree cannot be estimated using the PMSF model without the site frequencies file)

location = "dayhoff"
if (location == "local"){
  alignment_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/01_Data_all/"
  output_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/"
  repo_dir <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"
  iqtree2 <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_Software_IQ-Tree/iqtree-2.2.0-MacOSX/bin/iqtree2"
  
  number_parallel_processes <- 1
  
} else if (location == "soma"){
  alignment_dir <- "/data/caitlin/metazoan-mixtures/data_all/"
  output_dir <- "/data/caitlin/metazoan-mixtures/output/"
  repo_dir <- "/data/caitlin/metazoan-mixtures/"
  iqtree2 <- "/data/caitlin/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2"
  
  number_parallel_processes <- 2
  
} else if (location == "dayhoff"){
  alignment_dir <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/data_all/"
  output_dir <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/output/"
  repo_dir <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/"
  iqtree2 <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2"
  
  number_parallel_processes <- 3
  
} else if (location == "rona"){
  alignment_dir <- "/home/caitlin/metazoan-mixtures/data_all/"
  output_dir <- "/home/caitlin/metazoan-mixtures/output/"
  repo_dir <- "/home/caitlin/metazoan-mixtures/"
  iqtree2 <- "/home/caitlin/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2"
  
  number_parallel_processes <- 3
  
}

# Set parameters that are identical for all run locations
pmsf_initial_model <- "'LG+F+G'" # select guide tree from the C20, C60, LG+C20 and LG+C60 runs
pmsf_model <- c("'LG+C60+F+R4'", "'LG+C20+F+R4'", "'C60+F+R4'", "'C20+F+R4'") # extract best model from the C20, C60, LG+C20 and LG+C60 runs
pmsf_model_code <- c("PMSF_LG_C60", "PMSF_LG_C20", "PMSF_C60", "PMSF_C20")
iqtree_num_threads <- 15
ml_tree_bootstraps <- 1000
run_IQTREE <- FALSE
trim_incomplete <- TRUE



#### 2. Prepare functions, variables and packages ####
# Open packages
library(parallel)
library(ape)

# Source functions
source(paste0(repo_dir, "code/func_pmsf_trees.R"))

# Create a new folder to store the pmsf trees in
pmsf_dir <- paste0(output_dir, "pmsf_trees/")
if (dir.exists(pmsf_dir) == FALSE){dir.create(pmsf_dir)}

# Prepare file paths for tsv files
df_op_01_01 <- paste0(output_dir, "01_01_maximum_likelihood_tree_estimation_parameters.tsv")
df_op_pmsf_params <- paste0(output_dir, "01_05_PMSF_tree_estimation_parameters.tsv")
df_op_pmsf_commands <- paste0(output_dir, "01_06_PMSF_output.tsv")
pmsf_iqtree_command_files <- paste0(output_dir, "01_05_PMSF_iqtree_commands_", 1:3, ".txt")
pmsf_unrun_iqtree_command_files <- paste0(output_dir, "01_05_PMSF_iqtree_commands_", 1:3, "_toRun.txt")
pmsf_test_output_files <- paste0(output_dir, "01_06_PMSF_output_step",1:3,".tsv")

# Extract the list of all files from the folder containing alignments/models/partition files
all_files <- list.files(alignment_dir, recursive = TRUE)
# Remove any folder with "00" in the name (extra files associated with alignments)
all_files <- grep("00_", all_files, value = TRUE, invert = TRUE)
# Add full file path
if (length(all_files) > 0){
  all_files <- paste0(alignment_dir, all_files)
}

# Extract the list of alignments (i.e. files that contain the word "alignment")
all_alignments <- grep("\\.alignment\\.", all_files, value = T)



#### 3. Estimate a tree for each alignment using the PMSF model in IQ-Tree ####
if (file.exists(df_op_pmsf_params) == TRUE){
  pmsf_df <- read.table(file = df_op_pmsf_params, header = TRUE, sep = "\t")
} else if (file.exists(df_op_pmsf_params) == FALSE){
  ## Prepare the parameters to estimate the PMSF trees
  # Open the maximum likelihood tree estimation parameters tsv
  ml_df <- read.table(file = df_op_01_01, header = TRUE, stringsAsFactors = FALSE)
  # Reduce only to the MFP rows
  ml_df <- ml_df[ml_df$model_code == "ModelFinder", ]
  # Sort maximum likelihood df by alignment name
  ml_df <- ml_df[order(ml_df$dataset),]
  # Construct a grid of the alignments and the four model codes
  pmsf_df <- expand.grid(ml_df$dataset, pmsf_model_code)
  names(pmsf_df) <- c("dataset", "model_code")
  pmsf_df$dataset <- as.character(pmsf_df$dataset)
  pmsf_df$model_code <- as.character(pmsf_df$model_code)
  # Add actual models as column
  pmsf_df$guide_tree_model <- pmsf_initial_model
  pmsf_df$pmsf_model <- rep(pmsf_model, each = length(all_alignments))
  # Sort pmsf_df by dataset name
  pmsf_df <- pmsf_df[order(pmsf_df$dataset),]
  # Update columns for PMSF run
  pmsf_df$matrix_name <- rep(unique(ml_df$matrix_name), each = length(pmsf_model_code))
  pmsf_df$sequence_format <- rep(unique(ml_df$sequence_format), each = length(pmsf_model_code))
  pmsf_df$alignment_file <- rep(sort(all_alignments), each = length(pmsf_model_code))
  pmsf_df$prefix <- paste0(pmsf_df$dataset, ".", pmsf_df$matrix_name, ".", pmsf_df$model_code)
  pmsf_df$iqtree_num_bootstraps <- ml_tree_bootstraps
  pmsf_df$iqtree_num_threads <- iqtree_num_threads
  pmsf_df$iqtree_path <- iqtree2
  pmsf_df$pmsf_dir <- pmsf_dir
  # Reorder columns
  pmsf_df <- pmsf_df[, c("dataset", "model_code", "guide_tree_model", "pmsf_model", "matrix_name", "prefix", "sequence_format", "iqtree_num_threads", "iqtree_num_bootstraps", 
                         "alignment_file", "iqtree_path", "pmsf_dir")]
  # Redo row names
  row.names(pmsf_df) <- 1:nrow(pmsf_df)
  # Save dataframe as output
  write.table(pmsf_df, file = df_op_pmsf_params, row.names = FALSE, sep = "\t")
}

# Change to the PMSF tree directory
setwd(pmsf_dir)

# Exclude Simion 2017 from the runs
pmsf_df <- pmsf_df[pmsf_df$dataset != "Simion2017", ]
# Renumber rows 
rownames(pmsf_df) <- 1:nrow(pmsf_df)

# 1. Estimate guide tree under simple model
guide_list <- lapply(1:nrow(pmsf_df), estimate.guide.tree.wrapper, pmsf_parameter_dataframe = pmsf_df, run.iqtree = FALSE)
guide_df <- as.data.frame(do.call(rbind, guide_list))
names(guide_df) <- c("IQTree_command_1", "guide_tree_suffix")
pmsf_command_df <- cbind(pmsf_df, guide_df)
# Paste the pmsf directory to the guide tree suffix
pmsf_command_df$guide_tree_path <- paste0(pmsf_dir, pmsf_command_df$guide_tree_suffix, ".treefile")
# Write the IQ-Tree commands out to a file
write(unique(pmsf_command_df$IQTree_command_1), file = pmsf_iqtree_command_files[1])
# Write the dataframe out to a file (to capture the output at this step)
write.table(pmsf_command_df, file = pmsf_test_output_files[1], row.names = FALSE, sep = "\t")
# Get the unrun pmsf_commands for estimating guide trees
unrun_guide_trees <- unique(pmsf_command_df$IQTree_command_1[!file.exists(pmsf_command_df$guide_tree_path)])
# Write the unrun IQ-Tree commands out to a file
if (length(unrun_guide_trees) > 0){
  write(unrun_guide_trees, file = pmsf_unrun_iqtree_command_files[1])
  if (run_IQTREE == TRUE){
    mclapply(unrun_guide_trees, system, mc.cores = number_parallel_processes)
  }
}

# 2. Perform the first phase of the PMSF model: estimate mixture model parameters given the guide tree and infer site-specific 
#   frequency profile (printed to .sitefreq file)
sitefreq_list <- lapply(1:nrow(pmsf_command_df), output.site.frequency.file.wrapper, pmsf_parameter_dataframe = pmsf_command_df, run.iqtree = FALSE)
sitefreq_df <- as.data.frame(do.call(rbind, sitefreq_list))
names(sitefreq_df) <- c( "IQTree_command_2", "site_frequencies_suffix")
pmsf_command_df <- cbind(pmsf_command_df, sitefreq_df)
# Paste the pmsf directory to the site tree suffix
pmsf_command_df$site_frequencies_path <- paste0(pmsf_dir, pmsf_command_df$site_frequencies_suffix, ".sitefreq")
# Write the IQ-Tree commands out to a file
write(unique(pmsf_command_df$IQTree_command_2), file = pmsf_iqtree_command_files[2])
# Write the dataframe out to a file (to capture the output at this step)
write.table(pmsf_command_df, file = pmsf_test_output_files[2], row.names = FALSE, sep = "\t")
# Get the unrun pmsf_commands for estimating guide trees
unrun_ssfp <- unique(pmsf_command_df$IQTree_command_2[!file.exists(pmsf_command_df$site_frequencies_path)])
# Write the unrun IQ-Tree commands out to a file
if (length(unrun_ssfp) > 0){
  write(unrun_ssfp, file = pmsf_unrun_iqtree_command_files[2])
  if (run_IQTREE == TRUE){
    mclapply(unique(pmsf_command_df$IQTree_command_2), system, mc.cores = number_parallel_processes)
  }
}


# 3. Perform the second phase of the PMSF model: conduct typical analysis using the inferred frequency model (instead of the mixture model) 
#   to save RAM and running time.
if (trim_incomplete == TRUE){
  trimmed_pcdf <- pmsf_command_df[file.exists(pmsf_command_df$site_frequencies_path),]
  pmsf_command_df <- trimmed_pcdf
}
tree_list <- lapply(1:nrow(pmsf_command_df), estimate.tree.with.inferred.PMSF.model.wrapper, pmsf_parameter_dataframe = pmsf_command_df, run.iqtree = FALSE)
tree_df <- as.data.frame(do.call(rbind, tree_list))
names(tree_df) <- c("IQTree_command_3", "pmsf_prefix")
pmsf_command_df <- cbind(pmsf_command_df, tree_df)
# Write the IQ-Tree commands out to a file
write(unique(pmsf_command_df$IQTree_command_3), file = pmsf_iqtree_command_files[3])
# Write the dataframe out to a file (to capture the output at this step)
write.table(pmsf_command_df, file = pmsf_test_output_files[3], row.names = FALSE, sep = "\t")
if (run_IQTREE == TRUE){
  mclapply(unique(pmsf_command_df$IQTree_command_3), system, mc.cores = number_parallel_processes)
}

# 4. Extract the PMSF files for each run
op_files_list <- lapply(pmsf_command_df$pmsf_prefix, find.pmsf.files, pmsf_dir)
op_files_df <- as.data.frame(do.call(rbind, op_files_list))
names(op_files_df) <- c("pmsf_iqtree_file", "pmsf_log_file", "pmsf_tree_file")
pmsf_command_df <- cbind(pmsf_command_df, op_files_df)

# Save the dataframe of completed PMSF tree information
write.table(pmsf_command_df, file = df_op_pmsf_commands, row.names = FALSE, sep = "\t")


