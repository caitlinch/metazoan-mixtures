# metazoan-mixtures/code/01_estimate_PMSF_trees.R
## This script estimates maximum likelihood trees and maximum likelihood under constraint trees for empirical data sets
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
# alignment_dir       <- Directory containing alignments for all data sets
#                           Alignment naming convention: [manuscript].[matrix_name].[sequence_type].fa
#                           E.g. Cherryh2022.alignment1.aa.fa
# output_dir          <- Directory for IQ-Tree output (trees and tree mixtures)
# repo_dir            <- Location of caitlinch/metazoan-mixtures github repository
# iqtree2             <- Location of IQ-Tree2 stable release

# pmsf_initial_model        <- Model used to estimate guide tree for the site-specific frequency model (we use pmsf_initial_model = "'LG+C60+F+G'")
#                                 Place model inside two different types of quotes so it can be pasted into an IQ-Tree command line properly, e.g. "'model+R6'"
# iqtree_num_threads        <- Number of parallel threads for IQ-Tree to use. Can be a set number (e.g. 2) or "AUTO"
# ml_tree_bootstraps        <- Number of ultrafast bootstraps (UFB) to perform in IQ-Tree
# number_parallel_processes <- The number of simultaneous processes to run at once using mclapply(). 
#                                 If 1, then all processes will run sequentially

location = "local"
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
  
  number_parallel_processes <- 20
  
} else if (location == "dayhoff"){
  alignment_dir <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/data_all/"
  output_dir <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/output/"
  repo_dir <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/"
  iqtree2 <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2"
  
  number_parallel_processes <- 1
  
} else if (location == "rona"){
  alignment_dir <- "/home/caitlin/metazoan-mixtures/data_all/"
  output_dir <- "/home/caitlin/metazoan-mixtures/output/"
  repo_dir <- "/home/caitlin/metazoan-mixtures/"
  iqtree2 <- "/home/caitlin/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2"
  
  number_parallel_processes <- 1
  
}

# Set parameters that are identical for all run locations
pmsf_initial_model <- "'LG+C60+F+G'"
iqtree_num_threads <- 5
ml_tree_bootstraps <- 1000



#### 2. Prepare functions, variables and packages ####
# Open packages
library(parallel)
library(ape)

# Source functions
source(paste0(repo_dir, "code/func_pmsf_trees.R"))
# source(paste0(repo_dir, "code/func_estimate_trees.R"))
# source(paste0(repo_dir, "code/func_data_processing.R"))

# Create a new folder to store the pmsf trees in
pmsf_dir <- paste0(output_dir, "pmsf_trees/")
if (dir.exists(pmsf_dir) == FALSE){dir.create(pmsf_dir)}

# Prepare file paths for tsv files
df_op_01_01 <- paste0(output_dir, "01_01_maximum_likelihood_tree_estimation_parameters.tsv")
df_op_pmsf_params <- paste0(output_dir, "01_05_PMSF_tree_estimation_parameters.tsv")
df_op_pmsf_commands <- paste0(output_dir, "01_06_PMSF_iqtree_commands.tsv")

# Extract the list of all files from the folder containing alignments/models/partition files
all_files <- list.files(alignment_dir, recursive = TRUE)
if (length(all_files) > 0){
  all_files <- paste0(alignment_dir, all_files)
}

# Extract the list of alignments (i.e. files that contain the word "alignment")
all_alignments <- grep("\\.alignment\\.", all_files, value = T)

alignment_path <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/02_maximum_likelihood_trees/00_pmsf_tests/Philippe2011.UPDUNN_MB_FixedNames.aa.alignment.nex"
alignment_path <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/02_maximum_likelihood_trees/00_pmsf_tests/Dunn2008.Dunn2008_FixedNames.aa.alignment.fasta"
simple_model = "'LG+C20+F+G'"
alignment_prefix = "Dunn2008"
iqtree_path = iqtree2
pmsf_dir = "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/02_maximum_likelihood_trees/00_pmsf_tests/"



#### 3. Estimate a tree for each alignment using the PMSF model in IQ-Tree ####
## Prepare the parameters to estimate the PMSF trees
# Open the maximum likelihood tree estimation parameters tsv
pmsf_df <- read.table(file = df_op_01_01, header = TRUE, stringsAsFactors = FALSE)
# Reduce only to the MFP rows
pmsf_df <- pmsf_df[pmsf_df$model_code == "ModelFinder", ]
# Update columns for PMSF run
pmsf_df$model_code <- "PMSF"
pmsf_df$prefix <- paste0(pmsf_df$dataset, ".", pmsf_df$matrix_name, ".", pmsf_df$model_code)
pmsf_df$guide_tree_model <- pmsf_initial_model
pmsf_df$iqtree_num_bootstraps <- ml_tree_bootstraps
pmsf_df$iqtree_num_threads <- iqtree_num_threads
pmsf_df$iqtree_path <- iqtree2
pmsf_df$pmsf_dir <- pmsf_dir
# Remove unnecessary columns
pmsf_df <- pmsf_df[, c("dataset", "model_code", "guide_tree_model", "matrix_name", "prefix", "sequence_format", "iqtree_num_threads", "iqtree_num_bootstraps", 
                       "alignment_file", "iqtree_path", "pmsf_dir")]
# Save dataframe as output
write.table(pmsf_df, file = df_op_pmsf_params, row.names = FALSE, sep = "\t")

## Apply the PMSF wrapper function to determine the iqtree calls for tree estimation
pmsf_command_list <- lapply(1:nrow(pmsf_df), estimate.PMSF.tree.wrapper, pmsf_df, run.iqtree = FALSE)
pmsf_command_df <- as.data.frame(do.call(rbind, pmsf_command_list))
# Save dataframe as output
write.table(pmsf_command_df, file = df_op_pmsf_commands, row.names = FALSE, sep = "\t")



