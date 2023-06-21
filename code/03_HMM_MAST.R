## caitlinch/metazoan-mixtures/code/03_HMM_MAST.R
# This script applies the MAST (Mixtures Across Sites and Trees) method to constrained hypothesis trees estimated from 14 empirical Metazoan datasets
# Caitlin Cherryh 2023


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
  mast_dir                <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/05_all_MAST_output/"
  au_test_dir             <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/05_au_test/"
  output_dir              <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/01_output_files/"
  repo_dir                <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"
  iqtree2                 <- "iqtree2"
  iqtree_tm               <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_Software_IQ-Tree/iqtree-2.2.0.8.mix.1.hmm-MacOSX/bin/iqtree2"
  iqtree_hmmster          <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_Software_IQ-Tree/iqtree-2.2.3.hmmster-MacOSX/bin/iqtree"
  
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
  iqtree_hmmster          <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.3.hmmster-Linux/bin/iqtree2"
  
  ## Phylogenetic and IQ-Tree2 parameters
  iqtree_num_threads      <- 10
}

## Control parameters
run.MAST                  <- FALSE
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



#### 4. Apply mixtures across trees and sites (MAST model) - phyloHMM ####
# Create phyloHMM command lines in IQ-Tree
phyloHMM_run_list <- lapply(1:nrow(model_df), phyloHMM.wrapper, mast_df = model_df, iqtree_tree_mixtures = iqtree_tm,
                            MAST_branch_length_option = "TR", iqtree_num_threads = iqtree_num_threads, iqtree_min_branch_length = 0.00001,
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
if (run.MAST == TRUE){
  # Call IQ-Tree2
  system(phyloHMM_df$phyloHMM_iqtree2_command)
}



#### 5. Apply mixtures across trees and sites (MAST model) - HMMster ####
# Create HMMster command lines in IQ-Tree
HMMster_run_list <- lapply(1:nrow(model_df), HMMster.wrapper, mast_df = model_df, iqtree_tree_mixtures = iqtree_hmmster, 
                           MAST_branch_length_option = "T", iqtree_num_threads = iqtree_num_threads, iqtree_min_branch_length = 0.00001,
                           run.iqtree = FALSE)
HMMster_run_df <- as.data.frame(do.call(rbind, HMMster_run_list))
# Bind dataframe
HMMster_df <- cbind(model_df,HMMster_run_df)
# Write dataframe
HMMster_call_path <- paste0(output_dir, "03_02_HMMster_command_lines.tsv")
write.table(HMMster_df, file = HMMster_call_path, sep = "\t")
# Write command lines as text file
HMMster_call_text_path <- paste0(output_dir, "03_02_HMMster_command_lines.txt")
write(HMMster_df$HMMster_iqtree2_command, file = HMMster_call_text_path)
# Run HMMster
if (run.MAST == TRUE){
  # Call IQ-Tree2
  system(HMMster_df$HMMster_iqtree2_command)
}



#### 6. Extract results from phyloHMM and HMMster runs ####
## Extract file paths
# Extract list of HMMster and phyloHMM files
all_files <- list.files(mast_dir, recursive = TRUE)
phylohmm_files <- grep("\\.hmm", grep("phyloHMM", all_files, value = T), value = T)
hmmster_files <- grep("\\.hmm", grep("HMMster", all_files, value = T), value = T)
## HMM results and site ratios
# To extract information from the completed HMM run:
phylohmm_hmm_output <- lapply(paste0(mast_dir, phylohmm_files), extract.HMM.output)
# To extract information from the completed HMM run:
hmmster_hmm_output <- lapply(paste0(mast_dir, hmmster_files), extract.HMM.output)
## Tree Weights
# To extract information from the tree weights:
phylohmm_tw_output <- lapply(paste0(mast_dir, gsub("\\.hmm", ".iqtree", phylohmm_files)), extract.tree.weights)
# To extract information from the tree weights:
hmmster_tw_output <- lapply(paste0(mast_dir, gsub("\\.hmm", ".iqtree", hmmster_files)), extract.tree.weights)
## Create HMM dataframe
hmm_df <- rbind(as.data.frame(do.call(rbind, phylohmm_hmm_output)), as.data.frame(do.call(rbind, hmmster_hmm_output)))
## Create tree weight dataframe
tw_df <- rbind(as.data.frame(do.call(rbind, phylohmm_tw_output)), as.data.frame(do.call(rbind, hmmster_tw_output)))
## Format dataframes
hmm_df$dataset <- unlist(lapply(1:nrow(hmm_df), function(i){strsplit(hmm_df$hmm_file[i], "\\.")[[1]][1]}))
hmm_df$matrix_name <- unlist(lapply(1:nrow(hmm_df), function(i){strsplit(hmm_df$hmm_file[i], "\\.")[[1]][2]}))
hmm_df$model_code <- unlist(lapply(1:nrow(hmm_df), function(i){strsplit(hmm_df$hmm_file[i], "\\.")[[1]][3]}))
hmm_df <- hmm_df[, c("dataset", "matrix_name", "model_code",  "analysis_type", "number_hypothesis_trees",
                     "tree_1_hmm_probs", "tree_2_hmm_probs", "tree_3_hmm_probs", "tree_4_hmm_probs",
                     "tree_5_hmm_probs", "tree_1_number_sites", "tree_2_number_sites", "tree_3_number_sites",
                     "tree_4_number_sites", "tree_5_number_sites", "tree_1_ratio_sites", "tree_2_ratio_sites",
                     "tree_3_ratio_sites", "tree_4_ratio_sites", "tree_5_ratio_sites")]
names(hmm_df) <- c("hmm_dataset", "hmm_matrix_name", "hmm_model_code",  "hmm_analysis_type", "hmm_number_hypothesis_trees",
                   "tree_1_hmm_probs", "tree_2_hmm_probs", "tree_3_hmm_probs", "tree_4_hmm_probs",
                   "tree_5_hmm_probs", "tree_1_number_sites", "tree_2_number_sites", "tree_3_number_sites",
                   "tree_4_number_sites", "tree_5_number_sites", "tree_1_ratio_sites", "tree_2_ratio_sites",
                   "tree_3_ratio_sites", "tree_4_ratio_sites", "tree_5_ratio_sites")
hmm_df <- hmm_df[order(hmm_df$hmm_analysis_type, hmm_df$hmm_dataset, hmm_df$hmm_matrix_name, hmm_df$hmm_matrix_name), ]
tw_df$dataset <- unlist(lapply(1:nrow(tw_df), function(i){strsplit(tw_df$iq_file[i], "\\.")[[1]][1]}))
tw_df$matrix_name <- unlist(lapply(1:nrow(tw_df), function(i){strsplit(tw_df$iq_file[i], "\\.")[[1]][2]}))
tw_df$model_code <- unlist(lapply(1:nrow(tw_df), function(i){strsplit(tw_df$iq_file[i], "\\.")[[1]][3]}))
tw_df <- tw_df[, c("dataset", "matrix_name", "model_code",  "analysis_type", "number_hypothesis_trees",
                   "tree_1_tree_weight", "tree_2_tree_weight", "tree_3_tree_weight", "tree_4_tree_weight",
                   "tree_5_tree_weight", "tree_1_total_tree_length", "tree_2_total_tree_length",
                   "tree_3_total_tree_length", "tree_4_total_tree_length", "tree_5_total_tree_length",
                   "tree_1_sum_internal_branch_lengths", "tree_2_sum_internal_branch_lengths",
                   "tree_3_sum_internal_branch_lengths", "tree_4_sum_internal_branch_lengths",
                   "tree_5_sum_internal_branch_lengths")]
names(tw_df) <- c("tw_dataset", "tw_matrix_name", "tw_model_code",  "tw_analysis_type", "tw_number_hypothesis_trees",
                  "tree_1_tree_weight", "tree_2_tree_weight", "tree_3_tree_weight", "tree_4_tree_weight",
                  "tree_5_tree_weight", "tree_1_total_tree_length", "tree_2_total_tree_length",
                  "tree_3_total_tree_length", "tree_4_total_tree_length", "tree_5_total_tree_length",
                  "tree_1_sum_internal_branch_lengths", "tree_2_sum_internal_branch_lengths",
                  "tree_3_sum_internal_branch_lengths", "tree_4_sum_internal_branch_lengths",
                  "tree_5_sum_internal_branch_lengths")
tw_df <- tw_df[order(tw_df$tw_analysis_type, tw_df$tw_dataset, tw_df$tw_matrix_name, tw_df$tw_matrix_name), ]
## Combine dataframes
mast_df <- cbind(hmm_df, tw_df)
mast_df <- mast_df[, c("hmm_dataset", "hmm_matrix_name", "hmm_model_code",  "hmm_analysis_type", "hmm_number_hypothesis_trees",
                       "tree_1_hmm_probs", "tree_2_hmm_probs", "tree_3_hmm_probs", "tree_4_hmm_probs",
                       "tree_5_hmm_probs", "tree_1_number_sites", "tree_2_number_sites", "tree_3_number_sites",
                       "tree_4_number_sites", "tree_5_number_sites", "tree_1_ratio_sites", "tree_2_ratio_sites",
                       "tree_3_ratio_sites", "tree_4_ratio_sites", "tree_5_ratio_sites",
                       "tree_1_tree_weight", "tree_2_tree_weight", "tree_3_tree_weight", "tree_4_tree_weight",
                       "tree_5_tree_weight", "tree_1_total_tree_length", "tree_2_total_tree_length",
                       "tree_3_total_tree_length", "tree_4_total_tree_length", "tree_5_total_tree_length",
                       "tree_1_sum_internal_branch_lengths", "tree_2_sum_internal_branch_lengths",
                       "tree_3_sum_internal_branch_lengths", "tree_4_sum_internal_branch_lengths",
                       "tree_5_sum_internal_branch_lengths")]
names(mast_df) <- c("dataset", "matrix_name", "model_code",  "analysis_type", "number_hypothesis_trees",
                    "tree_1_hmm_probs", "tree_2_hmm_probs", "tree_3_hmm_probs", "tree_4_hmm_probs",
                    "tree_5_hmm_probs", "tree_1_number_sites", "tree_2_number_sites", "tree_3_number_sites",
                    "tree_4_number_sites", "tree_5_number_sites", "tree_1_ratio_sites", "tree_2_ratio_sites",
                    "tree_3_ratio_sites", "tree_4_ratio_sites", "tree_5_ratio_sites",
                    "tree_1_tree_weight", "tree_2_tree_weight", "tree_3_tree_weight", "tree_4_tree_weight",
                    "tree_5_tree_weight", "tree_1_total_tree_length", "tree_2_total_tree_length",
                    "tree_3_total_tree_length", "tree_4_total_tree_length", "tree_5_total_tree_length",
                    "tree_1_sum_internal_branch_lengths", "tree_2_sum_internal_branch_lengths",
                    "tree_3_sum_internal_branch_lengths", "tree_4_sum_internal_branch_lengths",
                    "tree_5_sum_internal_branch_lengths")
## Save output dataframe
mast_df_file <- paste0(output_dir, "04_01_HMM_MAST_model_output.tsv")
write.table(mast_df, mast_df_file, sep= "\t")



#### 7. Apply AU test to each dataset ####
# Run the tree topology tests
if (run.tree.topology.tests == TRUE){
  top_test_call_list <- lapply(1:nrow(phyloHMM_df), tree.topology.test.wrapper, df = phyloHMM_df, output_dir = au_test_dir, iqtree2 = iqtree2, iqtree_num_threads = iqtree_num_threads,
                               iqtree_num_RELL_replicates = 10000, run.iqtree = TRUE)
  au_test_calls <- unlist(top_test_call_list)
  write(au_test_calls, paste0(output_dir, "03_02_au_test_calls.text"))
}
# Extract the tree topology test results
all_op_files <- list.files(au_test_dir)
au_test_iqtree_files <- paste0(au_test_dir, grep("AU_test", grep("\\.iqtree", all_op_files, value = TRUE), value = TRUE))
au_test_list <- lapply(au_test_iqtree_files, extract.tree.topology.test.results)
# Save tree topology test results to file
au_test_df <- as.data.frame(do.call(rbind, au_test_list))
write.table(au_test_df, paste0(output_dir, "04_01_tree_topology_test_results.tsv"), sep= "\t")


