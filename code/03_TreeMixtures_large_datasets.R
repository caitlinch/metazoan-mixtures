## caitlinch/metazoan-mixtures/code/04_TreeMixtures_large_datasets.R
# This script applies the MAST model to two large datasets (Simion 2017 and Hejnol 2009)
# Caitlin Cherryh 2023


#### 1. Input parameters ####
## Specify parameters:
# alignment_dir               <- Directory containing alignments for all data sets
#                                   Alignment naming convention: [manuscript].[matrix_name].[sequence_type].fa
#                                   E.g. Cherryh2022.alignment1.aa.fa\
# big_data_output_dir         <- Directory for constraint trees, hypothesis trees, and MAST run for the two large datasets
# output_dir                  <- Directory for IQ-Tree output (trees and tree mixtures)
# repo_dir                    <- Location of caitlinch/metazoan-mixtures github repository

# iqtree2                     <- Location of IQ-Tree2 stable release
# iqtree_tm                   <- Location of IQ-Tree2 phyloHMM release (currently: version 2.2.0.8.mix.1.hmm)
# iqtree_hmmster              <- Location of IQ-Tree2 HMMster release (currently: version 2.2.3.hmmster)
# iqtree_num_threads          <- Number of parallel threads for IQ-Tree to use. Can be a set number (e.g. 2) or "AUTO"
#                                     See IQ-Tree documentation for more details (http://www.iqtree.org/doc/Command-Reference)
# hypothesis_tree_bootstraps  <- Number of ultrafast bootstraps (UFB) to perform in IQ-Tree when estimating constrained maximum likelihood trees
# number_parallel_processes   <- The number of simultaneous processes to run at once using mclapply(). 
#                                     If 1, then all processes will run sequentially

## Specify control parameter values (all take logical values TRUE or FALSE):
# prepare.hypothesis.trees      <- TRUE to prepare constraint trees and create command lines to estimate hypothesis trees (constrained maximum likelihood trees). FALSE to skip.
# estimate.hypothesis.trees     <- TRUE to estimate all hypothesis trees (constrained maximum likelihood trees). FALSE to skip.
# run.MAST.model                <- TRUE to call IQ-Tree2 and run the MAST model. FALSE to output IQ-Tree2 command lines without running MAST model.
# run.tree.topology.tests       <- TRUE to call IQ-Tree2 and run the tree topology tests. FALSE to output IQ-Tree2 command lines without running tree topology tests.

location = "local"
if (location == "local"){
  alignment_dir         <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/01_Data_all/"
  big_data_output_dir   <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/05_large_datasets/"
  output_dir            <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/01_output_files/"
  repo_dir              <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"
  
  iqtree2               <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_Software_IQ-Tree/iqtree-2.2.0-MacOSX/bin/iqtree2"
  iqtree_tm             <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_Software_IQ-Tree/iqtree-2.2.0.8.mix.1.hmm-MacOSX/bin/iqtree2"
  iqtree_hmmster        <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_Software_IQ-Tree/iqtree-2.2.3.hmmster-MacOSX/bin/iqtree"
  iqtree_num_threads        <- 3
  hypothesis_tree_bootstraps <- 1000
  number_parallel_processes <- 1
} else if (location == "dayhoff"){
  alignment_dir         <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/data_all/"
  big_data_output_dir   <- ""
  output_dir            <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/output/"
  repo_dir              <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/"
  
  iqtree2               <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2"
  iqtree_tm               <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0.8.mix.1.hmm-Linux/bin/iqtree2"
  iqtree_hmmster          <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.3.hmmster-Linux/bin/iqtree2"
  iqtree_num_threads        <- 20
  hypothesis_tree_bootstraps <- 1000
  number_parallel_processes <- 4
} 

# Set control parameters
control_parameters <- list(prepare.hypothesis.trees = FALSE,
                           estimate.hypothesis.trees = FALSE,
                           run.MAST.model = FALSE,
                           run.tree.topology.tests = FALSE)



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

# Create file paths for output files
output_file_paths <- paste0(output_dir, c())



#### 3. Prepare constraint trees ####



#### 4. Construct constraint trees ####



#### 5. Estimate hypothesis trees ####



#### 6. Apply MAST model ####



#### 7. Apply tree topology tests ####




