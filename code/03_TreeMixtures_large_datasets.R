## caitlinch/metazoan-mixtures/code/04_TreeMixtures_large_datasets.R
# This script applies the MAST model to two large datasets (Simion 2017 and Hejnol 2009)
# Caitlin Cherryh 2023

# The analyses in this script refer to the following papers:
## Hejnol 2009
# Paper:
# Hejnol, A., Obst, M., Stamatakis, A. 2009, "Assessing the root of bilaterian animals with scalable phylogenomic methods", 
# Proc. R. Soc. B., 276:4261-4270, http://doi.org/10.1098/rspb.2009.0896
# Dataset:
# Available as supplementary information at http://doi.org/10.1098/rspb.2009.0896


## Simion 2017
# Paper: 
# Simion, P., Philippe, H., Baurain, D., et al 2017, "A Large and Consistent Phylogenomic Dataset Supports Sponges as the Sister Group to All Other Animals",
# Current Biology, 27:958-967, https://doi.org/10.1016/j.cub.2017.02.031.
# Dataset:
# Simion, P. 2017, "SuppData_Metazoa_2017", Available on GitHub at https://github.com/psimion/SuppData_Metazoa_2017

# Analysis notes:
#   - Skip ultrafast bootstraps for these large datasets (time constraint)


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
# number_parallel_processes   <- The number of simultaneous processes to run at once using mclapply(). 
#                                     If 1, then all processes will run sequentially

## Specify control parameter values (all take logical values TRUE or FALSE):
# prepare.hypothesis.trees      <- TRUE to prepare constraint trees and create command lines to estimate hypothesis trees (constrained maximum likelihood trees). FALSE to skip.
# estimate.hypothesis.trees     <- TRUE to estimate all hypothesis trees (constrained maximum likelihood trees). FALSE to skip.
# run.MAST.model                <- TRUE to call IQ-Tree2 and run the MAST model. FALSE to output IQ-Tree2 command lines without running MAST model.
# run.tree.topology.tests       <- TRUE to call IQ-Tree2 and run the tree topology tests. FALSE to output IQ-Tree2 command lines without running tree topology tests.

location = "soma"
if (location == "local"){
  alignment_dir         <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/01_Data_all/"
  big_data_output_dir   <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/05_large_datasets/"
  output_dir            <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/01_output_files/"
  repo_dir              <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"
  
  iqtree2               <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_Software_IQ-Tree/iqtree-2.2.0-MacOSX/bin/iqtree2"
  iqtree_tm             <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_Software_IQ-Tree/iqtree-2.2.0.8.mix.1.hmm-MacOSX/bin/iqtree2"
  iqtree_hmmster        <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_Software_IQ-Tree/iqtree-2.2.3.hmmster-MacOSX/bin/iqtree"
  iqtree_num_threads        <- 3
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
  number_parallel_processes <- 4
} else if (location == "soma"){
  alignment_dir         <- "/home/caitlin/metazoan-mixtures/data_all/"
  big_data_output_dir   <- "/home/caitlin/metazoan-mixtures/large_datasets/"
  output_dir            <- "/home/caitlin/metazoan-mixtures/output_csvs/"
  repo_dir              <- "/home/caitlin/metazoan-mixtures/"
  
  iqtree2               <- "/home/caitlin/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2"
  iqtree_tm               <- "/home/caitlin/metazoan-mixtures/iqtree/iqtree-2.2.0.8.mix.1.hmm-Linux/bin/iqtree2"
  iqtree_hmmster          <- "/home/caitlin/metazoan-mixtures/iqtree/iqtree-2.2.3.hmmster-Linux/bin/iqtree2"
  iqtree_num_threads        <- 20
  number_parallel_processes <- 4
}



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
output_file_paths <- c(paste0(output_dir, c("01_02_maximum_likelihood_included_taxa.tsv", "LargeDatasets_hypothesis_tree_estimation_calls.txt")),
                       paste0(big_data_output_dir, c("Simion2017.supermatrix_97sp_401632pos_1719genes.constrained_ML.hypothesis_trees.treefile",
                                                     "Hejnol2009.Hejnol_etal_2009_FixedNames.constrained_ML.hypothesis_trees.treefile")) )



#### 3. Construct constraint trees ####
alignment_taxa_df <- read.table(output_file_paths[[1]], header = T)
simion_constraint_trees <- output.constraint.trees(dataset = "Simion2017", matrix_name = "supermatrix_97sp_401632pos_1719genes", 
                                                   constraint_tree_dir = big_data_output_dir, dataset_info = all_datasets,
                                                   matrix_taxa_info = matrix_taxa, ml_tree_tips_df = alignment_taxa_df, 
                                                   force.update.constraint.trees = TRUE)
hejnol_constraint_trees <- output.constraint.trees(dataset = "Hejnol2009", matrix_name = "Hejnol_etal_2009_FixedNames", 
                                                   constraint_tree_dir = big_data_output_dir, dataset_info = all_datasets,
                                                   matrix_taxa_info = matrix_taxa, ml_tree_tips_df = alignment_taxa_df, 
                                                   force.update.constraint.trees = TRUE)



#### 4. Estimate hypothesis trees ####
## Models of substitution from original papers
# Simion dataset: partitioned by gene, each gene has independent LG+G4+F model
# Hejnol dataset: "Models of molecular evolution were evaluated using the Perl script available from the RAxML website. 
#                  ML searches and bootstrap analyses were executed under the Gamma model of rate heterogeneity" (Hejnol et al. 2009)

### Construct model search for Hejnol 2009 partitions
# Prepare file paths and iqtree parameters
hejnol_al_file <- paste0(alignment_dir, grep("Hejnol2009", list.files(alignment_dir), value = T))
hejnol_partition_file <- paste0(big_data_output_dir, grep("gene_partitions", grep("Hejnol2009", list.files(big_data_output_dir), value = T), value = T))
hejnol_partition_prefix <- paste0(paste(strsplit(basename(hejnol_partition_file), "\\.")[[1]][1:2], collapse = "."), ".models")
# Assemble iqtree command line
hejnol_partition_model_call <- paste0(iqtree2, " -s ", hejnol_al_file, " -spp ", hejnol_partition_file, " -m TESTMERGEONLY ", 
                                      " -nt ", iqtree_num_threads, " -pre ", hejnol_partition_prefix)

### Estimate hypothesis trees
## Simion 2017 dataset
# Prepare file paths and iqtree parameters
simion_al_file <- paste0(alignment_dir, grep("Simion2017", list.files(alignment_dir), value = T))
simion_partition_file <- paste0(big_data_output_dir, grep("gene_partitions.models", grep("Simion2017", list.files(big_data_output_dir), value = T), value = T))
simion_constraint_trees <- sort(paste0(big_data_output_dir, grep("constraint_tree", grep("Simion2017", list.files(big_data_output_dir), value = T), value = T)))
simion_hypothesis_prefixes <- paste0("Simion2017.supermatrix_97sp_401632pos_1719genes.ML_H", 1:length(simion_constraint_trees))
# Assemble iqtree command line
# # With partitions - do not use this for MAST runs
# simion_hypothesis_tree_calls <- paste0(iqtree2, " -s ", simion_al_file, " -spp ", simion_partition_file, " -g ", simion_constraint_trees, 
#                                        " -nt ", iqtree_num_threads, " -pre ", simion_hypothesis_prefixes)
# Without partitions - concatenated LG+G4+F model
simion_hypothesis_tree_calls <- paste0(iqtree2, " -s ", simion_al_file, " -m LG+G4+F -g ", simion_constraint_trees, 
                                       " -nt ", iqtree_num_threads, " -pre ", simion_hypothesis_prefixes)
## Hejnol 2009 dataset
# Prepare file paths and iqtree parameters
hejnol_al_file <- paste0(alignment_dir, grep("Hejnol2009", list.files(alignment_dir), value = T))
hejnol_partition_file <- paste0(big_data_output_dir, grep("gene_partitions.models", grep("Hejnol2009", list.files(big_data_output_dir), value = T), value = T))
hejnol_constraint_trees <- sort(paste0(big_data_output_dir, grep("constraint_tree", grep("Hejnol2009", list.files(big_data_output_dir), value = T), value = T)))
hejnol_hypothesis_prefixes <- paste0("Hejnol2009.Hejnol_etal_2009_FixedNames.ML_H", 1:length(hejnol_constraint_trees))
# Assemble iqtree command line
hejnol_hypothesis_tree_calls <- paste0(iqtree2, " -s ", hejnol_al_file, " -spp ", hejnol_partition_file, " -g ", hejnol_constraint_trees, 
                                       " -nt ", iqtree_num_threads, " -pre ", hejnol_hypothesis_prefixes)

## Output all command lines
write(c(simion_hypothesis_tree_calls, hejnol_hypothesis_tree_calls), file = output_file_paths[[2]])

## Collate trees
# Collect all files from directory
all_large_dataset_files <- list.files(big_data_output_dir)
# Collate Simion 2017 trees into single file
simion_trees <- sort(paste0(big_data_output_dir, grep("ML_H", grep("Simion2017", all_large_dataset_files, value = T), value = T)))
simion_collated_trees <- c(unlist(lapply(simion_trees, readLines)), "")
write(simion_collated_trees, file = output_file_paths[[3]])
# Collate Hejnol 2009 trees into single file
hejnol_trees <- sort(paste0(big_data_output_dir, grep("ML_H", grep("Hejnol2009", all_large_dataset_files, value = T), value = T)))
hejnol_collated_trees <- c(unlist(lapply(hejnol_trees, readLines)), "")
write(hejnol_collated_trees, file = output_file_paths[[4]])


#### 5. Apply MAST model ####
# NOTE: Before running MAST model for the Hejnol dataset, the model of substitution for the MAST run must be set. This model cannot be a partition model!

# Filepaths for iqtree commands
simion_al_file <- paste0(alignment_dir, grep("Simion2017", list.files(alignment_dir), value = T))
simion_partition_file <- paste0(big_data_output_dir, grep("gene_partitions.models", grep("Simion2017", list.files(big_data_output_dir), value = T), value = T))
simion_collated_trees <- output_file_paths[[3]]
hejnol_al_file <- paste0(alignment_dir, grep("Hejnol2009", list.files(alignment_dir), value = T))
hejnol_partition_file <- paste0(big_data_output_dir, grep("gene_partitions.models", grep("Hejnol2009", list.files(big_data_output_dir), value = T), value = T))
hejnol_collated_trees <- output_file_paths[[4]]
# Parameters to set minimum branch lengths:
simion_num_sites = 401632
hejnol_num_sites = 270580

#### phyloHMM
#     $ iqtree2 -m "TMIX{GTR+FO+G,GTR+FO+G}+TR" -hmm -te data1.all.top.txt -s data1.fa -blmin 0.00001 -nt 30 -pre phylohmm
## Simion 2017 dataset
simion_phylohmm_call <- paste0(iqtree_tm, " -m 'LG+G4+F+TR' -hmm -te ", simion_collated_trees, " -s ", simion_al_file, 
                               " -blmin ", format(1/simion_num_sites, scientific = F, digits = 1),
                              " -nt ", iqtree_num_threads, " -pre Simion2017.supermatrix_97sp_401632pos_1719genes.phyloHMM")
## Hejnol 2009 dataset
hejnol_phylohmm_call <- paste0(iqtree_tm, " -m 'LG+G4+F+TR' -hmm -te ", hejnol_collated_trees, " -s ", hejnol_al_file, 
                               " -blmin ", format(1/hejnol_num_sites, scientific = F, digits = 1),
                               " -nt ", iqtree_num_threads, " -pre Hejnol2009.Hejnol_etal_2009_FixedNames.phyloHMM")

#### HMMster
#     $ iqtree2 -m "TMIX{GTR+FO+G,GTR+FO+G}+T" -te hypothesis_trees.treefile -s alignment.fa  -hmmster -blmin 0.00001 -nt 30 -pre hmmster
## Simion 2017 dataset
simion_hmmster_call <- paste0(iqtree_hmmster, " -m 'LG+G4+F+T' -te ", simion_collated_trees, " -s ", simion_al_file, 
                              " -hmmster -blmin ", format(1/simion_num_sites, scientific = F, digits = 1),
                              " -nt ", iqtree_num_threads, " -pre Simion2017.supermatrix_97sp_401632pos_1719genes.HMMster")
## Hejnol 2009 dataset
# Assemble iqtree command line
hejnol_hmmster_call <- paste0(iqtree_hmmster, " -m 'LG+G4+F+T' -te ", hejnol_collated_trees, " -s ", hejnol_al_file, 
                              " -hmmster -blmin ", format(1/hejnol_num_sites, scientific = F, digits = 1),
                              " -nt ", iqtree_num_threads, " -pre Hejnol2009.Hejnol_etal_2009_FixedNames.HMMster")



#### 6. Apply tree topology tests ####
## Simion 2017 dataset
simion_au_test_call <- paste0(iqtree2, " -s ", simion_al_file, " -spp ", simion_partition_file, " -n 0 -z ", simion_collated_trees, " -zb 1000 -au -zw",
                              " -nt ", iqtree_num_threads, " -pre Simion2017.supermatrix_97sp_401632pos_1719genes.AU_test")
## Hejnol 2009 dataset
hejnol_au_test_call <- paste0(iqtree2, " -s ", hejnol_al_file, " -spp ", hejnol_partition_file, " -n 0 -z ", hejnol_collated_trees, " -zb 1000 -au -zw",
                              " -nt ", iqtree_num_threads, " -pre Hejnol2009.Hejnol_etal_2009_FixedNames.AU_test")



