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

# iqtree_num_threads        <- Number of parallel threads for IQ-Tree to use. Can be a set number (e.g. 2) or "AUTO"
# iqtree_mrate              <- Specify a comma separated list of rate heterogeneity types for model selection in IQ-Tree
#                                 We set iqtree_mrate = "E,I,G,I+G,R,I+R"
#                                 See IQ-Tree documentation for more details (http://www.iqtree.org/doc/Command-Reference)
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
  
}

# Set parameters that are identical for all run locations
iqtree_mrate <- "E,I,G,I+G,R,I+R"
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

# Extract the list of all files from the folder containing alignments/models/partition files
all_files <- list.files(alignment_dir)
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







