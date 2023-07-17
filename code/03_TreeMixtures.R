## caitlinch/metazoan-mixtures/code/03_TreeMixtures.R
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
# mast_dir                <- Directory for MAST output
# output_dir              <- Directory for output csvs
# repo_dir                <- Location of caitlinch/metazoan-mixtures github repository
# iqtree2                 <- Location of IQ-Tree2 stable release (version 2.2.2)
# iqtree_MAST             <- Location of IQ-Tree2 MAST release (currently: version 2.2.3.hmmster)

## Phylogenetic and IQ-Tree2 parameters
# iqtree_num_threads      <- Number of parallel threads for IQ-Tree to use. Can be a set number (e.g. 2) or "AUTO"

## Control parameters
# run.MAST                  <- TRUE to call IQ-Tree2 and run the MAST model. FALSE to output IQ-Tree2 command lines without running MAST model.
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
  iqtree_MAST             <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_Software_IQ-Tree/iqtree-2.2.3.hmmster-MacOSX/bin/iqtree"
  
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
  iqtree_MAST             <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.3.hmmster-Linux/bin/iqtree2"
  
  ## Phylogenetic and IQ-Tree2 parameters
  iqtree_num_threads      <- 25
} else if (location == "rona"){
  ## File paths
  alignment_dir           <- "/home/caitlin/metazoan-mixtures/data_all/"
  hypothesis_tree_dir     <- "/home/caitlin/metazoan-mixtures/hypothesis_trees/"
  pmsf_sitefreq_dir       <- "/home/caitlin/metazoan-mixtures/pmsf_sitefreqs/"
  mast_dir                <- "/home/caitlin/metazoan-mixtures/mast/"
  au_test_dir             <- "/home/caitlin/metazoan-mixtures/au_test/"
  output_dir              <- "/home/caitlin/metazoan-mixtures/output_csvs/"
  repo_dir                <- "/home/caitlin/metazoan-mixtures/"
  iqtree2                 <- "/home/caitlin/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2"
  iqtree_MAST             <- "/home/caitlin/metazoan-mixtures/iqtree/iqtree-2.2.3.hmmster-Linux/bin/iqtree2"
  
  ## Phylogenetic and IQ-Tree2 parameters
  iqtree_num_threads      <- 50
}

## Control parameters
control_parameters <- list(prepare.MAST = FALSE,
                           run.MAST = FALSE,
                           extract.MAST = FALSE,
                           prepare.tree.topology.tests = FALSE,
                           run.tree.topology.tests = FALSE,
                           extract.tree.topology.tests = FALSE)



#### 2. Prepare variables, open packages and source functions ####
# Extend the number of digits allowed (so BIC and logL can be extracted properly from iqtree files)
options(digits = 12)

# Source files containing functions
source(paste0(repo_dir, "code/func_data_processing.R"))
source(paste0(repo_dir, "code/func_estimate_trees.R"))



#### 3. Prepare for analysis ####
## Check whether the dataframe for phyloHMM runs has been created
mast_parameter_path <- paste0(output_dir, "03_01_MAST_parameters.tsv")
if (file.exists(mast_parameter_path) == TRUE){
  model_df <- read.table(mast_parameter_path, header = T)
} else if (file.exists(mast_parameter_path) == FALSE){
  ## Prepare the dataframe 
  model_df <- read.table(paste0(output_dir, "01_03_best_models_per_alignment.tsv"), header = T)
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
  
  ## Remove all models without complete hypothesis trees
  model_df <- model_df[which(model_df$model_class == "profile"), ]
  
  ## Prepare the hypothesis tree files
  model_df$hypothesis_tree_path <- unlist(lapply(paste0(model_df$dataset, ".", model_df$matrix_name, ".", model_df$model_code), 
                                                 collate.hypothesis.trees, input_dir = hypothesis_tree_dir,
                                                 output_dir = mast_dir))
  ## Add the  file paths to the dataframe
  model_df$best_model_sitefreq_path <- paste0(pmsf_sitefreq_dir, basename(model_df$best_model_sitefreq_path))
  model_df$alignment_path <- paste0(alignment_dir, basename(model_df$alignment_path))
  model_df$hypothesis_tree_path <- paste0(hypothesis_tree_dir, basename(model_df$hypothesis_tree_path))
  
  ## Add column with the minimum branch length for each alignment
  al_dims_file <- paste0(output_dir, grep("summary_alignment_details", list.files(output_dir), value = T))
  if (file.exists(al_dims_file) == TRUE){
    # Open alignment dimensions file
    al_dims_df <- read.csv(al_dims_file, stringsAsFactors = F)
    # Reorder row order to match model_df
    include_rows <- match(model_df$matrix_name, al_dims_df$matrix_name)[which(!is.na(match(model_df$matrix_name, al_dims_df$matrix_name)))]
    al_dims_df <- al_dims_df[include_rows, ]
    # Sort dataframes by column value
    al_dims_df <- al_dims_df[order(al_dims_df$dataset, al_dims_df$matrix_name),]
    model_df <- model_df[order(model_df$dataset, model_df$matrix_name),]
    # Determine minimum branch length as 1/N for each dataset (where N is the total number of sites), and add it to the model_df as a column
    model_df$num_sites <- al_dims_df$num_sites
    model_df$min_MAST_bl_from_alignment <- round(100/model_df$num_sites, 5)
    model_df$min_MAST_bl_arbitrary <- 0.00001
  }
  
  ## Sort and remove columns
  model_df <- model_df[, c("dataset", "model_code", "matrix_name", "alignment_path", "sequence_format", "prefix", "best_model",                
                           "best_model_sitefreq_path", "best_model_LogL", "best_model_BIC", "best_model_wBIC", "tree_LogL",
                           "tree_UnconstrainedLogL", "tree_NumFreeParams", "tree_BIC", "tree_length", "tree_SumInternalBranch",
                           "tree_PercentInternalBranch",  "estimated_rates", "estimated_gamma", "estimated_state_frequencies",
                           "hypothesis_tree_path", "num_sites", "min_MAST_bl_from_alignment", "min_MAST_bl_arbitrary")]
  
  ## Write the dataframe
  write.table(model_df, file = mast_parameter_path, sep = "\t")
}




#### 4. Apply mixtures across trees and sites (MAST model) ####
if (control_parameters$prepare.MAST == TRUE){
  # Update parameter file paths for MAST run on server
  model_df$alignment_path <- paste0(alignment_dir, basename(model_df$alignment_path))
  model_df$best_model_sitefreq_path <- paste0(pmsf_sitefreq_dir, basename(model_df$best_model_sitefreq_path))
  model_df$hypothesis_tree_path <- paste0(hypothesis_tree_dir, basename(model_df$hypothesis_tree_path))
  
  # Create MAST command lines in IQ-Tree
  MAST_TR_run_df <- as.data.frame(do.call(rbind, lapply(1:nrow(model_df), function(i){MAST.wrapper(i, mast_df = model_df, 
                                                                                                   iqtree_tree_mixtures = iqtree_MAST, 
                                                                                                   MAST_branch_length_option = "TR",
                                                                                                   iqtree_num_threads = iqtree_num_threads,
                                                                                                   run.iqtree = FALSE) }) ) )
  MAST_T_run_df <- as.data.frame(do.call(rbind, lapply(1:nrow(model_df), function(i){MAST.wrapper(i, mast_df = model_df, 
                                                                                                  iqtree_tree_mixtures = iqtree_MAST, 
                                                                                                  MAST_branch_length_option = "T",
                                                                                                  iqtree_num_threads = iqtree_num_threads,
                                                                                                  run.iqtree = FALSE) }) ) )
  # Add branch length option column
  MAST_TR_run_df$MAST_branch_length_model <- "TR"
  MAST_T_run_df$MAST_branch_length_model <- "T"
  # Bind dataframe
  MAST_run_df <- rbind(cbind(model_df, MAST_TR_run_df), cbind(model_df, MAST_T_run_df))
  MAST_run_df$prefix <- paste0(MAST_run_df$prefix, ".", MAST_run_df$MAST_branch_length_model)
  # Write dataframe
  MAST_call_path <- paste0(output_dir, "03_02_MAST_command_lines.tsv")
  write.table(MAST_run_df, file = MAST_call_path, sep = "\t")
  # Write command lines as text file
  MAST_call_text_path <- paste0(output_dir, "03_02_MAST_command_lines.txt")
  write(MAST_run_df$MAST_iqtree2_command, file = MAST_call_text_path)
  # Run HMMster
  if (control_parameters$run.MAST == TRUE){
    # Call IQ-Tree2
    system(MAST_run_df$MAST_iqtree2_command)
  }
}


#### 5. Extract results from MAST runs ####
if (control_parameters$extract.MAST == TRUE){
  ## Extract file paths
  # Extract list of MAST files
  all_files <- list.files(mast_dir, recursive = TRUE)
  exclude_dirs <- "minbl_0.00001|HMM|HMMster|phyloHMM"
  op_files <- all_files[grep(exclude_dirs, dirname(all_files), invert = T)]
  mast_files <- grep("\\.iqtree", grep("MAST", op_files, value = T), value = T)
  
  ## Tree Weights
  # To extract information from the tree weights:
  mast_tws_df <- as.data.frame(do.call(rbind, lapply(mast_files, extract.tree.weights)))
  
  ## Format dataframes
  mast_tws_df$dataset <- unlist(lapply(1:nrow(mast_tws_df), function(i){strsplit(mast_tws_df$iq_file[i], "\\.")[[1]][1]}))
  mast_tws_df$matrix_name <- unlist(lapply(1:nrow(mast_tws_df), function(i){strsplit(mast_tws_df$iq_file[i], "\\.")[[1]][2]}))
  mast_tws_df$model_code <- unlist(lapply(1:nrow(mast_tws_df), function(i){strsplit(mast_tws_df$iq_file[i], "\\.")[[1]][3]}))
  mast_tws_df$mast_branch_type <- unlist(lapply(hmmster_files, return.MAST.branch.type))
  mast_tws_df$minimum_branch_length <- paste0("0.", unlist(lapply(1:nrow(mast_tws_df), function(i){strsplit(mast_tws_df$iq_file[i], "\\.")[[1]][7]})))
  mast_tws_df <- mast_tws_df[, c("dataset", "matrix_name", "model_code",  "analysis_type", 
                                 "mast_branch_type", "minimum_branch_length", "number_hypothesis_trees",
                                 "tree_1_tree_weight", "tree_2_tree_weight", "tree_3_tree_weight", "tree_4_tree_weight",
                                 "tree_5_tree_weight", "tree_1_total_tree_length", "tree_2_total_tree_length",
                                 "tree_3_total_tree_length", "tree_4_total_tree_length", "tree_5_total_tree_length",
                                 "tree_1_sum_internal_branch_lengths", "tree_2_sum_internal_branch_lengths",
                                 "tree_3_sum_internal_branch_lengths", "tree_4_sum_internal_branch_lengths",
                                 "tree_5_sum_internal_branch_lengths")]
  names(mast_tws_df) <- c("dataset", "matrix_name", "model_code",  "analysis_type", 
                          "mast_branch_type", "minimum_branch_length", "number_hypothesis_trees",
                          "tree_1_tree_weight", "tree_2_tree_weight", "tree_3_tree_weight", "tree_4_tree_weight",
                          "tree_5_tree_weight", "tree_1_total_tree_length", "tree_2_total_tree_length",
                          "tree_3_total_tree_length", "tree_4_total_tree_length", "tree_5_total_tree_length",
                          "tree_1_sum_internal_branch_lengths", "tree_2_sum_internal_branch_lengths",
                          "tree_3_sum_internal_branch_lengths", "tree_4_sum_internal_branch_lengths",
                          "tree_5_sum_internal_branch_lengths")
  mast_tws_df <- mast_tws_df[order(mast_tws_df$analysis_type, mast_tws_df$mast_branch_type, 
                                   mast_tws_df$dataset, mast_tws_df$matrix_name, mast_tws_df$model_code), ]
  
  ## Save output dataframe
  mast_df_file <- paste0(output_dir, "04_01_MAST_model_output.tsv")
  write.table(mast_tws_df, mast_df_file, sep= "\t")
}



#### 7. Apply AU test to each dataset ####
if (control_parameters$prepare.tree.topology.tests == TRUE){
  # Run the tree topology tests
  if (control_parameters$run.tree.topology.tests == TRUE){
    top_test_call_list <- lapply(1:nrow(model_df), tree.topology.test.wrapper, df = model_df, output_dir = au_test_dir, 
                                 iqtree2 = iqtree2, iqtree_num_threads = iqtree_num_threads,
                                 iqtree_num_RELL_replicates = 10000, run.iqtree = TRUE)
    au_test_calls <- unlist(top_test_call_list)
    write(au_test_calls, paste0(output_dir, "03_02_au_test_calls.text"))
  }
}

if (control_parameters$extract.tree.topology.tests == TRUE){
  # Extract the tree topology test results
  all_op_files <- list.files(au_test_dir)
  au_test_iqtree_files <- paste0(au_test_dir, grep("AU_test", grep("\\.iqtree", all_op_files, value = TRUE), value = TRUE))
  au_test_list <- lapply(au_test_iqtree_files, extract.tree.topology.test.results)
  # Save tree topology test results to file
  au_test_df <- as.data.frame(do.call(rbind, au_test_list))
  write.table(au_test_df, paste0(output_dir, "04_01_tree_topology_test_results.tsv"), sep= "\t")
}



