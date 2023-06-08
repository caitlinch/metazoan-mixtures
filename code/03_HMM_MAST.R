# metazoan-mixtures/code/02_HMM_MAST.R
## This script applies the MAST (Mixtures Across Sites and Trees) method to constrained hypothesis trees estimated from 14 empirical Metazoan datasets
# Caitlin Cherryh, 2023

# In this script, MAST refers to the Mixtures Across Sites and Trees model
#   Thomas KF Wong, Caitlin Cherryh, Allen G Rodrigo, Matthew W Hahn, Bui Quang Minh, Robert Lanfear 2022, 
#   "MAST: Phylogenetic Inference with Mixtures Across Sites and Trees", bioRxiv 2022.10.06.511210; 
#   doi: https://doi.org/10.1101/2022.10.06.511210



#### 1. Input parameters ####
## File paths
# alignment_dir           <- Directory containing alignments for all data sets
#                               Alignment naming convention: [manuscript].[matrix_name].[sequence_type].fa
#                               E.g. Cherryh2022.alignment1.aa.fa
# hypothesis_tree_dir     <- Directory containing constrained maximum likelihood trees
# pmsf_sitefreq_dir       <- Directory containing output files (including .sitefreq files) from IQ-Tree runs with PMSF model
# mast_dir                <- Directory for phyloHMM/MAST output
# output_dir              <- Directory for output csvs
# repo_dir                <- Location of caitlinch/metazoan-mixtures github repository
# iqtree2                 <- Location of IQ-Tree2 stable release (version 2.2.2)
# iqtree_tm               <- Location of IQ-Tree2 MAST release (currently: version 2.2.0.8.mix.1.hmm)

## Phylogenetic and IQ-Tree2 parameters
# iqtree_num_threads      <- Number of parallel threads for IQ-Tree to use. Can be a set number (e.g. 2) or "AUTO"

## Control parameters
# run.phyloHMM              <- TRUE to call IQ-Tree2 and run the phyloHMM. FALSE to output IQ-Tree2 command lines without running phyloHMM.
# run.tree.topology.tests   <- TRUE to call IQ-Tree2 and run the tree topology tests. FALSE to output IQ-Tree2 command lines without running tree topology tests.

location = "local"
if (location == "local"){
  ## File paths
  alignment_dir           <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/01_Data_all/"
  hypothesis_tree_dir     <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/04_hypothesis_trees/hypothesis_tree_estimation/"
  pmsf_sitefreq_dir       <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/02_pmsf_site_freqs/"
  mast_dir                <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/06_phyloHMM/"
  au_test_dir             <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/05_au_test/"
  output_dir              <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/"
  repo_dir                <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"
  iqtree2                 <- "iqtree2"
  iqtree_tm               <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_Software_IQ-Tree/iqtree-2.2.0.8.mix.1.hmm-MacOSX/bin/iqtree2"
  
  ## Phylogenetic and IQ-Tree2 parameters
  iqtree_num_threads      <- 3
} else if (location == "dayhoff"){
  ## File paths
  alignment_dir           <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/data_all/"
  hypothesis_tree_dir     <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/output/hyp_tree_output_files/"
  pmsf_sitefreq_dir       <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/output/pmsf_trees/"
  mast_dir                <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/output/phyloHMM/"
  au_test_dir             <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/output/au_test/"
  output_dir              <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/output/"
  repo_dir                <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/"
  iqtree2                 <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2"
  iqtree_tm               <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0.8.mix.1.hmm-Linux/bin/iqtree2"
  
  ## Phylogenetic and IQ-Tree2 parameters
  iqtree_num_threads      <- 10
}

## Control parameters
run.phyloHMM              <- FALSE
run.tree.topology.tests   <- FALSE



#### 2. Prepare variables, open packages and source functions ####
# Extend the number of digits allowed (so BIC and logL can be extracted properly from iqtree files)
options(digits = 12)

# Source files containing functions
source(paste0(repo_dir, "code/func_data_processing.R"))
source(paste0(repo_dir, "code/func_estimate_trees.R"))



#### 3. Prepare for analysis ####
## Check whether the dataframe for phyloHMM runs has been created
phylohmm_parameter_path <- paste0(output_dir, "03_01_phyloHMM_parameters.tsv")
if (file.exists(phylohmm_parameter_path) == TRUE){
  model_df <- read.table(phylohmm_parameter_path, header = T)
} else if (file.exists(phylohmm_parameter_path) == FALSE){
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
  ## Prepare the hypothesis tree files
  model_df$hypothesis_tree_path <- unlist(lapply(paste0(model_df$dataset, ".", model_df$matrix_name, ".", model_df$model_code), 
                                                 collate.hypothesis.trees, input_dir = hypothesis_tree_dir,
                                                 output_dir = mast_dir))
  ## Add the  file paths to the dataframe
  model_df$best_model_sitefreq_path <- paste0(pmsf_sitefreq_dir, basename(model_df$best_model_sitefreq_path))
  model_df$alignment_path <- paste0(alignment_dir, basename(model_df$alignment_path))
  model_df$hypothesis_tree_path <- paste0(hypothesis_tree_dir, basename(model_df$hypothesis_tree_path))
  ## Write the dataframe
  write.table(model_df, file = phylohmm_parameter_path, sep = "\t")
}



#### 4. Apply mixtures across trees and sites (MAST model) ####
# Create phyloHMM command lines in IQ-Tree
phyloHMM_run_list <- lapply(1:nrow(model_df), phyloHMM.wrapper, mast_df = model_df, MAST_branch_length_option = "TR",
                            iqtree_tree_mixtures = iqtree_tm, iqtree_num_threads = iqtree_num_threads, iqtree_min_branch_length = 0.00001,
                            run.iqtree = FALSE)
phyloHMM_run_df <- as.data.frame(do.call(rbind, phyloHMM_run_list))
# Bind dataframe
phyloHMM_df <- cbind(model_df,phyloHMM_run_df)
# Write dataframe
phylohmm_call_path <- paste0(output_dir, "03_02_phyloHMM_command_lines.tsv")
write.table(phyloHMM_df, file = phylohmm_call_path, sep = "\t")
# Write command lines as text file
phylohmm_call_text_path <- paste0(output_dir, "03_02_phyloHMM_command_lines.txt")
write(phyloHMM_df$phyloHMM_iqtree2_command, file = phylohmm_call_text_path)
# Run phyloHMM
if (run.phyloHMM == TRUE){
  # Call IQ-Tree2
  system(phyloHMM_df$phyloHMM_iqtree2_command)
}

# To extract information from the completed HMM run:
hmm_output <- extract.phyloHMM.output(output_prefix = output_prefix, output_directory = "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/constraint_trees/00_test_phyloHMM/")



#### 5. Apply AU test to each dataset ####
# Run the tree topology tests
if (run.tree.topology.tests == TRUE){
  top_test_call_list <- lapply(1:nrow(phyloHMM_df), tree.topology.test.wrapper, df = phyloHMM_df, output_dir = au_test_dir, iqtree2 = iqtree2, iqtree_num_threads = iqtree_num_threads,
                               iqtree_num_RELL_replicates = 10000, run.iqtree = TRUE)
  au_test_calls <- unlist(top_test_call_list)
  write(au_test_calls, paste0(output_dir, "03_02_au_test_calls.text"))
}
# Extract the tree topology test results
top_test_list <- tree.topology.test.wrapper(1:nrow(phyloHMM_df), df = phyloHMM_df, output_dir = NA, iqtree2 = iqtree2, iqtree_num_threads = iqtree_num_threads,
                                            iqtree_num_RELL_replicates = 10000, run.iqtree = FALSE, return.AU.output = TRUE)
top_test_df <- as.data.frame(do.call(rbind, top_test_list))
# Write the tree topology results out
tree_top_path <- paste0(output_dir, "03_03_tree_topology_test_results.tsv")
write.table(top_test_df, file = tree_top_path, sep = "\t")


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


