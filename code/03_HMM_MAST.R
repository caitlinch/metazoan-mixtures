# metazoan-mixtures/code/02_HMM_MAST.R
## This script applies the MAST (Mixtures Across Sites and Trees) method to constrained hypothesis trees estimated from 14 empirical Metazoan datasets
# Caitlin Cherryh, 2022

# In this script, MAST refers to the Mixtures Across Sites and Trees model
#   Thomas KF Wong, Caitlin Cherryh, Allen G Rodrigo, Matthew W Hahn, Bui Quang Minh, Robert Lanfear 2022, 
#   "MAST: Phylogenetic Inference with Mixtures Across Sites and Trees", bioRxiv 2022.10.06.511210; 
#   doi: https://doi.org/10.1101/2022.10.06.511210



#### 1. Input parameters ####
## Specify parameters:
# alignment_dir           <- Directory containing alignments for all data sets
#                               Alignment naming convention: [manuscript].[matrix_name].[sequence_type].fa
#                               E.g. Cherryh2022.alignment1.aa.fa
# hypothesis_tree_dir     <- Directory containing constrained maximum likelihood trees
# pmsf_sitefreq_dir       <- Directory containing output files (including .sitefreq files) from IQ-Tree runs with PMSF model
# mast_dir                <- Directory for phyloHMM/MAST output
# repo_dir                <- Location of caitlinch/metazoan-mixtures github repository

# iqtree2                 <- Location of IQ-Tree2 stable release (version 2.2.2)
# iqtree_tm               <- Location of IQ-Tree2 MAST release (currently: version 2.2.0.8.mix.1.hmm)
# iqtree_num_threads      <- Number of parallel threads for IQ-Tree to use. Can be a set number (e.g. 2) or "AUTO"

location = "local"
if (location == "local"){
  alignment_dir           <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/01_Data_all/"
  hypothesis_tree_dir     <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/04_hypothesis_trees/"
  pmsf_sitefreq_dir       <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/02_pmsf_site_freqs/"
  mast_dir                <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/05_tree_mixtures/"
  repo_dir                <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"
  
  iqtree2                 <- "iqtree2"
  iqtree_tm               <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_Software_IQ-Tree/iqtree-2.2.0.8.mix.1.hmm-MacOSX/bin/iqtree2"
  iqtree_num_threads      <- 3
}



#### 2. Prepare variables, open packages and source functions ####




#### 3. Prepare for analysis ####
## Prepare the dataframe 
model_df <- read.table(paste0(repo_dir, "output/01_03_best_models_per_alignment.tsv"), header = T)
# Add the alignments to the dataframe
all_aldir_files <- list.files(alignment_dir)
al_files <- grep("alignment", grep("00_", all_aldir_files, value = T, invert = T), value = T)
model_df$alignment_path <- unlist(lapply(1:nrow(model_df), function(i){
  grep(model_df$dataset[i], grep(model_df$matrix_name[i], al_files, value = T), value = T)
  }))
# Split the best_model column into two columns: one for the best model to input into IQ-Tree, 
#     and one for the path to the .sitefreq file (present only if best model is a PMSF model)
split_best_model_str                <- strsplit(model_df$best_model, ":")
model_df$best_model                 <- unlist(lapply(1:nrow(model_df), function(x){split_best_model_str[[x]][1]}))
model_df$best_model_sitefreq_path   <- unlist(lapply(1:nrow(model_df), function(x){split_best_model_str[[x]][2]}))
# Reorder the columns
model_df <- model_df[,c("dataset", "model_code", "matrix_name", "alignment_path", "sequence_format", "prefix",
                        "best_model", "best_model_sitefreq_path", "best_model_LogL", "best_model_BIC",
                        "best_model_wBIC", "tree_LogL", "tree_UnconstrainedLogL", "tree_NumFreeParams",
                        "tree_BIC", "tree_length", "tree_SumInternalBranch", "tree_PercentInternalBranch",
                        "estimated_rates", "estimated_gamma", "estimated_state_frequencies",
                        "maximum_likelihood_tree")]
# Add the  file paths to the dataframe
model_df$best_model_sitefreq_path <- paste0(pmsf_sitefreq_dir, model_df$best_model_sitefreq_path)
model_df$alignment_path <- paste0(alignment_dir, model_df$alignment_path)

# Prepare the hypothesis tree files




#### 4. Apply mixtures across trees and sites (MAST model) ####
# ## To run the MAST model and determine the weights of each tree (deprecated - now using phyloHMM):
# # Create the tree mixture prefix and command lines
# ml_tree_df$MAST_prefix <- paste0(ml_tree_df$prefix,".TR")
# ml_tree_df$MAST_call <- lapply(1:nrow(ml_tree_df), tree.mixture.wrapper, iqtree_tm_path = iqtree2_tm,
#                                iqtree_num_threads = iqtree_num_threads, df = ml_tree_df)
# # Run the mixture of trees models
# mclapply(ml_tree_df$MAST_call, system, mc.cores = number_parallel_processes)

## To run the phyloHMM
## Note: possible to apply HMM model with PMSF model?
test_trees_file <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/constraint_trees/00_test_phyloHMM/Test.Nosenko2013.nonribosomal_9187_smatrix.LG.ML_Hypothesis_trees.treefile"
tree_file <- test_trees_file
alignment_file <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/01_Data_all/Nosenko2013.nonribosomal_9187_smatrix.aa.alignment.phy"
output_prefix <- "test.Nosenko2013.nonribo.LG.HMM"
model <- "LG"
MAST_model <- paste0(model, "+TR") # branch-length restricted MAST model: where a branch occurs in multiple treesm it is constrained to have the same length in each tree

# To run phyloHMM for the toy example
iqtree_hmm_command <- run.phyloHMM(tree_file = tree_file, alignment_file = alignment_file, MAST_model = MAST_model, output_prefix = output_prefix, 
                                   iqtree_phyloHMM = iqtree_tm, iqtree_num_threads = "AUTO", run.iqtree = FALSE)
# To extract information from the completed HMM run:
hmm_output <- extract.phyloHMM.output(output_prefix = output_prefix, output_directory = "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/constraint_trees/00_test_phyloHMM/")



#### 5. Apply AU test to each dataset ####




#### 6. Extract results from IQ-Tree output files ####
# Incomplete code: still need to extract results from MAST model and save to tsv file
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


