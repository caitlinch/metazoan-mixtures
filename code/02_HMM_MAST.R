# metazoan-mixtures/code/02_HMM_MAST.R
## This script applies the HMM MAST model to empirical phylogenetic datasets
# Caitlin Cherryh, 2022

## This script:
# 1. Applies the MAST (Mixtures Across Sites and Trees) method to constrained hypothesis trees estimated from empirical Metazoan datasets

# In this script, MAST refers to the Mixtures Across Sites and Trees model
#   Thomas KF Wong, Caitlin Cherryh, Allen G Rodrigo, Matthew W Hahn, Bui Quang Minh, Robert Lanfear 2022, 
#   "MAST: Phylogenetic Inference with Mixtures Across Sites and Trees", bioRxiv 2022.10.06.511210; 
#   doi: https://doi.org/10.1101/2022.10.06.511210



#### 1. Input parameters ####
## Specify parameters:
# alignment_dir       <- Directory containing alignments for all data sets
#                       Alignment naming convention: [manuscript].[matrix_name].[sequence_type].fa
#                       E.g. Cherryh2022.alignment1.aa.fa
# output_dir          <- Directory for IQ-Tree output (trees and tree mixtures)
# repo_dir            <- Location of caitlinch/metazoan-mixtures github repository

# iqtree2             <- Location of IQ-Tree2 stable release (version 2.2.2)
# iqtree_tm           <- Location of IQ-Tree2 MAST release (currently: version 2.2.0.8.mix.1.hmm)

# iqtree_num_threads  <- Number of parallel threads for IQ-Tree to use. Can be a set number (e.g. 2) or "AUTO"

location = "local"
if (location == "local"){
alignment_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/01_Data_all/"
output_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/"
repo_dir <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"

iqtree2 <- "iqtree2"
iqtree_tm <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_Software_IQ-Tree/iqtree-2.2.0.8.mix.1.hmm-MacOSX/bin/iqtree2"

iqtree_num_threads = 3
}



#### 2. Prepare variables, open packages and source functions ####
test_trees_file <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/constraint_trees/00_test_phyloHMM/Test.Nosenko2013.nonribosomal_9187_smatrix.LG.ML_Hypothesis_trees.treefile"
tree_file <- test_trees_file
alignment_file <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/01_Data_all/Nosenko2013.nonribosomal_9187_smatrix.aa.alignment.phy"
output_prefix <- "test.Nosenko2013.nonribo.LG.HMM"
model <- "LG"
MAST_model <- paste0(model, "+TR") # branch-length restricted MAST model: where a branch occurs in multiple treesm it is constrained to have the same length in each tree
 
#### 3. Apply mixtures across trees and sites (MAST model) ####
# Create a folder for the ml trees and move to that folder
m_tree_dir <- paste0(output_dir, "tree_mixtures/")
if (file.exists(m_tree_dir) == FALSE){dir.create(m_tree_dir)}
setwd(m_tree_dir)

# Create the tree mixture prefix and command lines
ml_tree_df$MAST_prefix <- paste0(ml_tree_df$prefix,".TR")
ml_tree_df$MAST_call <- lapply(1:nrow(ml_tree_df), tree.mixture.wrapper, iqtree_tm_path = iqtree2_tm,
                               iqtree_num_threads = iqtree_num_threads, df = ml_tree_df)

# Run the mixture of trees models
mclapply(ml_tree_df$MAST_call, system, mc.cores = number_parallel_processes)


############### Incomplete code: still need to extract results from MAST model and save to tsv file
# # Identify iqtree files from tree mixture run
# all_files <- paste0(a_tm_op_dir, list.files(a_tm_op_dir))
# all_iqtree_files <- grep("\\.iqtree", all_files, value = TRUE)
# tree_mixture_tr_iqfile <- grep(paste0(a_m_prefix,".TR"), all_iqtree_files, value = TRUE)
# 
# # Extract information about each tree mixture model run
# tr_results <- extract.tree.mixture.results(tree_mixture_file = tree_mixture_tr_iqfile, 
#                                            dataset = a_dataset, prefix = paste0(a_m_prefix,".TR"), 
#                                            model = m, best_model = best_model, tree_branch_option = "TR")
# 
# # Output results dataframe
# op_file <- paste0(a_tm_op_dir, a_m_prefix, "_tree_mixture_results.tsv")
# write.table(tr_results, file = op_file, row.names = FALSE, sep = "\t")

