## caitlinch/metazoan-mixtures/code/05_reformat_output_dataframes.R
# This script reads in csv/tsv files created prior in the pipeline, and reformats them nicely for the manuscript
# Caitlin Cherryh 2023


#### 1. Input parameters ####
## File paths
# output_file_dir         <- Directory for output csvs
# repo_dir                <- Location of caitlinch/metazoan-mixtures github repository
# hypothesis_tree_dir     <- Directory containing all constrained ML trees (i.e., the hypothesis trees)

## File paths
output_file_dir         <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/01_output_files/"
repo_dir                <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"
hypothesis_tree_dir     <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/04_hypothesis_trees/"



#### 2. Prepare variables, open packages and source functions ####
# Change the default number of digits (so we don't lose the decimal points in the BIC scores)
options(digits = 12)

# Open packages
library(readxl)

# Source files containing functions
source(paste0(repo_dir, "code/func_data_processing.R"))

# List all files in output directory
all_output_files <- paste0(output_file_dir, list.files(output_file_dir))
# # Remove any files with all 5 trees - only want to look at the output for the 2 tree model
# all_output_files <- grep("5trees", all_output_files, value = TRUE, invert = TRUE)



#### 3. Prepare summary of the manual topology check csv file ####
### Summarise topology results (as percentage of each output topology)
summary_topology_file <- paste0(output_file_dir, "summary_ML_tree_topology.csv")
if (file.exists(summary_topology_file) == FALSE){
  # Read in .xlsx file with manual topology check results
  topology_check_x_file <- grep("results_02_ML_tree_topology_ManualCheck", all_output_files, value = TRUE)
  topology_check_x_file <- grep("Summary", topology_check_x_file, value = TRUE, invert = TRUE)
  topology_check_df <- as.data.frame(read_excel(path = topology_check_x_file, sheet = "Topology"))
  # Remove Simion and Hejnol datasets - too computationally intensive to run full ML models
  topology_check_df <- topology_check_df[which(topology_check_df$dataset != "Hejnol2009" & topology_check_df$dataset != "Simion2017"), ]
  # Convert "NA" to NA
  topology_check_df$model_BIC[which(topology_check_df$model_BIC == "NA")] <- NA
  topology_check_df$sister_group[which(topology_check_df$sister_group == "NA")] <- NA
  topology_check_df$UFB_support[which(topology_check_df$UFB_support == "NA")] <- NA
  topology_check_df$PORI_topology[which(topology_check_df$PORI_topology == "NA")] <- NA
  topology_check_df$`CTEN+CNID_monophyletic`[which(topology_check_df$`CTEN+CNID_monophyletic` == "NA")] <- NA
  topology_check_df$PLAC_present[which(topology_check_df$PLAC_present == "NA")] <- NA
  topology_check_df$tree_BIC[which(topology_check_df$tree_BIC == "NA")] <- NA
  # Convert necessary columns to numeric
  topology_check_df$model_BIC <- as.numeric(topology_check_df$model_BIC)
  topology_check_df$tree_BIC <- as.numeric(topology_check_df$tree_BIC) 
  # Remove GTR20 model (had to use GTR_no_I for MAST, as GTR+I didn't run properly)
  topology_noGTR20_df <- topology_check_df[which(topology_check_df$model_code != "GTR20"), ]
  # Summarise results for each dataset as a percentage
  dataset_ids <- unique(paste0(topology_check_df$dataset, ".", topology_check_df$matrix_name))
  summary_topology_list <- lapply(dataset_ids, summarise.topology.results, topology_noGTR20_df, 
                                  excluded_models = c("C10", "C30", "C40", "C50"))
  summary_topology_df <-  as.data.frame(do.call(rbind, summary_topology_list))
  # Sort output by year
  summary_topology_df <- summary_topology_df[order(summary_topology_df$dataset_year, summary_topology_df$dataset, summary_topology_df$matrix_name),]
  # Write the output
  write.csv(summary_topology_df, file = summary_topology_file, row.names = FALSE)
}

### Nicely format output data frames of tree topologies and of sponge topologies
model_order <- c("PMSF_C20", "PMSF_C60", "PMSF_LG_C20", "PMSF_LG_C60", 
                 "C20", "C60", "LG_C20", "LG_C60", "CF4", "EHO", "EX_EHO",
                 "EX2", "EX3", "GTR20", "GTR_no_I", "JTT", "JTTDCMut", "LG", "LG4M",
                 "mtZOA", "PMB", "Poisson", "rtREV", "UL2", "UL3", "WAG",
                 "ModelFinder")
## For tree topologies:
tree_topology_file <- paste0(output_file_dir, "all_models_ML_tree_topology.csv")
if (file.exists(tree_topology_file) == FALSE){
  # Read in summary dataframe
  summary_topology_df <- read.csv(summary_topology_file, stringsAsFactors = FALSE)
  ordered_dataset_ids <- paste0(summary_topology_df$dataset, ".", summary_topology_df$matrix_name)
  # Extract tree topologies from dataframe
  tree_topology_list <- lapply(ordered_dataset_ids, tree.topology.results, topology_check_df, model_order)
  tree_topology_df <- as.data.frame(do.call(cbind, tree_topology_list))
  # Format dataframe
  names(tree_topology_df) <- ordered_dataset_ids
  tree_topology_df$row_names <- c("dataset", "matrix_name", model_order)
  tree_topology_df <- tree_topology_df[, c("row_names", ordered_dataset_ids)]
  # Output dataframe
  write.csv(tree_topology_df, file = tree_topology_file, row.names = FALSE)
}
## For Porifera topologies:
pori_topology_file <- paste0(output_file_dir, "all_models_ML_Porifera_topology.csv")
if (file.exists(pori_topology_file) == FALSE){
  # Read in summary dataframe
  summary_topology_df <- read.csv(summary_topology_file, stringsAsFactors = FALSE)
  ordered_dataset_ids <- paste0(summary_topology_df$dataset, ".", summary_topology_df$matrix_name)
  # Extract tree topologies from dataframe
  pori_topology_list <- lapply(ordered_dataset_ids, porifera.topology.results, topology_check_df, model_order)
  pori_topology_df <- as.data.frame(do.call(cbind, pori_topology_list))
  # Format dataframe
  names(pori_topology_df) <- ordered_dataset_ids
  pori_topology_df$row_names <- c("dataset", "matrix_name", model_order)
  pori_topology_df <- pori_topology_df[, c("row_names", ordered_dataset_ids)]
  # Output dataframe
  write.csv(pori_topology_df, file = pori_topology_file, row.names = FALSE)
}



#### 4. Prepare summary of the AU test results using the tree topology test tsv file ####
# Read in csv file
tree_test_df <- read.csv(file = grep("tree_topology_test_results.csv", all_output_files, value = TRUE), header = TRUE)
# Process for 2trees and 5trees analysis
analyses <- c("2_trees", "5_trees")
for (a in analyses){
  # Process each dataset one at a time
  temp_tree_test <- tree_test_df[tree_test_df$hypothesis_tree_analysis == a, ]
  summary_au_test_df <- as.data.frame(do.call(rbind, lapply(unique(temp_tree_test$ID), summarise.AU.test.results, temp_tree_test)))
  # Add new columns for model class and hypothesis tree analysis
  summary_au_test_df$model_class <- factor(summary_au_test_df$best_model_code,
                                           levels = c("LG_C60", "C60", "LG_C20", "PMSF_C60", "PMSF_LG_C60", "LG4M", "UL3"),
                                           labels = c("CXX", "CXX", "CXX", "PMSF", "PMSF", "Other", "Other"),
                                           ordered = TRUE)
  summary_au_test_df$hypothesis_tree_analysis <- a
  # Sort output by year
  summary_au_test_df <- summary_au_test_df[order(summary_au_test_df$year, summary_au_test_df$dataset, summary_au_test_df$matrix),]
  # Reorder columns
  summary_au_test_df <- summary_au_test_df[, c("hypothesis_tree_analysis", "dataset", "matrix", "model_class", "best_model_code",
                                               "topology_test", "tree_1", "tree_2", "tree_3", "tree_4", "tree_5", "year")]
  # Remove columns consisting only of NA
  summary_au_test_df <- Filter(function(x)!all(is.na(x)), summary_au_test_df)
  # Write the output
  summary_au_test_file <- paste0(output_file_dir, "summary_au_test_results_", gsub("_", "", a), ".csv")
  write.csv(summary_au_test_df, file = summary_au_test_file, row.names = FALSE)
}



#### 5. Prepare summary of the expected likelihood weights using the tree topology test tsv file ####
# Read in csv file
tree_test_df <- read.csv(file = grep("tree_topology_test_results.csv", all_output_files, value = TRUE), header = TRUE)
# Process for 2trees and 5trees analysis
analyses <- c("2_trees", "5_trees")
for (a in analyses){
  # Process each dataset one at a time
  temp_tree_test <- tree_test_df[tree_test_df$hypothesis_tree_analysis == a, ]
  # Process each dataset one at a time
  summary_elw_test_df <- as.data.frame(do.call(rbind, lapply(unique(temp_tree_test$ID), summarise.eLW, temp_tree_test)))
  # Add new column for model class
  summary_elw_test_df$model_class <- factor(summary_elw_test_df$best_model_code,
                                            levels = c("LG_C60", "C60", "LG_C20", "PMSF_C60", "PMSF_LG_C60", "LG4M", "UL3"),
                                            labels = c("CXX", "CXX", "CXX", "PMSF", "PMSF", "Other", "Other"),
                                            ordered = TRUE)
  summary_elw_test_df$hypothesis_tree_analysis <- a
  # Sort output by year
  summary_elw_test_df <- summary_elw_test_df[order(summary_elw_test_df$year, summary_elw_test_df$dataset, summary_elw_test_df$matrix),]
  # Reorder columns
  summary_elw_test_df <- summary_elw_test_df[, c("hypothesis_tree_analysis", "dataset", "matrix", "model_class", "best_model_code",
                                                 "topology_test", "tree_1", "tree_2", "tree_3", "tree_4", "tree_5", "year")]
  # Remove columns consisting only of NA
  summary_elw_test_df <- Filter(function(x)!all(is.na(x)), summary_elw_test_df)
  # Write the output
  summary_elw_test_file <- paste0(output_file_dir, "summary_elw_results_", gsub("_", "", a), ".csv")
  write.csv(summary_elw_test_df, file = summary_elw_test_file, row.names = FALSE)
}



#### 5. Prepare summary of the MAST tsv file ####
### Output site ratios from HMM weights
# Read in tsv file
mast_df <- read.csv(file = grep("MAST_model_output.csv", all_output_files, value = TRUE), header = TRUE)
# Process for 2trees and 5trees analysis
analyses <- c("2_trees", "5_trees")
for (a in analyses){
  # Output tree weights
  # Process each dataset one at a time
  temp_tree_test <- mast_df[mast_df$hypothesis_tree_analysis == a, ]
  mast_summary_df <- as.data.frame(do.call(rbind, lapply(1:nrow(temp_tree_test), summarise.tree.weights, tw_df = temp_tree_test)))
  # Sort output by year
  mast_summary_df <- mast_summary_df[order(mast_summary_df$year, mast_summary_df$dataset, mast_summary_df$matrix_name), ]
  mast_summary_df$hypothesis_tree_analysis <- a
  # Reorganise columns
  mast_summary_df <- mast_summary_df[,c("hypothesis_tree_analysis", "dataset", "matrix_name", "model_class", "model_code", 
                                        "tree_1_tree_weight", "tree_2_tree_weight", "tree_3_tree_weight", 
                                        "tree_4_tree_weight", "tree_5_tree_weight",
                                        "mast_branch_type", "minimum_branch_length", "number_hypothesis_trees", "year")]
  # Remove columns consisting only of NA
  mast_summary_df <- Filter(function(x)!all(is.na(x)), mast_summary_df)
  # Write the output for MAST tree weights
  write.csv(mast_summary_df, file = paste0(output_file_dir, "summary_MAST_treeWeight_results_", gsub("_", "", a), ".csv"), row.names = FALSE)
}



#### 6. Combine BIC from MAST and single tree ####
# Open MAST parameter and MAST output paths
mast_output_path <- grep("results_06_MAST_output", all_output_files, value = T)
mast_output <- read.csv(mast_output_path, stringsAsFactors = F)
ml_results_file <- grep("results_01_maximum_likelihood_iqtreeOutput", all_output_files, value = T)
ml_results <- read.csv(ml_results_file, stringsAsFactors = F)
ml_results$model_class <- factor(ml_results$model_class,
                                 levels = c("Single", "Other", "PMSF", "CXX"),
                                 labels = c("Q", "Mixture", "PMSF", "PM"),
                                 ordered = T)
# # To fix minbl lengths for MAST (if required_)
# mast_output$minimum_branch_length <- as.numeric(paste0("0.", unlist(lapply(strsplit(basename(mast_output$iq_file), "\\."), function(x){x[[grep("minbl", x)+1]]}))))
# Create parameters dataframe
datasets_df <- unique(mast_output[, c("dataset", "matrix_name")])
num_datasets = nrow(datasets_df)
datasets_df <- rbind(datasets_df, datasets_df, datasets_df, datasets_df)
datasets_df$model_class <- c(rep("Q", num_datasets), rep("Mixture", num_datasets), rep("PMSF", num_datasets), rep("PM", num_datasets))
rownames(datasets_df) <- 1:nrow(datasets_df)
# For each row in the datasets_df:
bic_list <- lapply(1:nrow(datasets_df), compare.multitree.BIC.wrapper, datasets_df = datasets_df, ml_results = ml_results, mast_output = mast_output)
bic_df <- as.data.frame(do.call(rbind, bic_list))
# Save BIC results
write.csv(bic_df, file = paste0(output_file_dir, "summary_BIC_values.csv"), row.names = FALSE)



#### 7. Combine log likelihood from MAST and single tree ####
# Open MAST parameter and MAST output paths
mast_output_path <- grep("results_06_MAST_output", all_output_files, value = T)
mast_output <- read.csv(mast_output_path, stringsAsFactors = F)
ml_results_file <- grep("results_01_maximum_likelihood_iqtreeOutput", all_output_files, value = T)
ml_results <- read.csv(ml_results_file, stringsAsFactors = F)
ml_results$model_class <- factor(ml_results$model_class,
                                 levels = c("Single", "Other", "PMSF", "CXX"),
                                 labels = c("Q", "Mixture", "PMSF", "PM"),
                                 ordered = T)
# # To fix minbl lengths for MAST (if required_)
# mast_output$minimum_branch_length <- as.numeric(paste0("0.", unlist(lapply(strsplit(basename(mast_output$iq_file), "\\."), function(x){x[[grep("minbl", x)+1]]}))))
# Create parameters dataframe
datasets_df <- unique(mast_output[, c("dataset", "matrix_name")])
num_datasets = nrow(datasets_df)
datasets_df <- rbind(datasets_df, datasets_df, datasets_df, datasets_df)
datasets_df$model_class <- c(rep("Q", num_datasets), rep("Mixture", num_datasets), rep("PMSF", num_datasets), rep("PM", num_datasets))
rownames(datasets_df) <- 1:nrow(datasets_df)
# For each row in the datasets_df:
logl_list <- lapply(1:nrow(datasets_df), compare.multitree.log.likelihood.wrapper, datasets_df = datasets_df, ml_results = ml_results, mast_output = mast_output)
logl_df <- as.data.frame(do.call(rbind, logl_list))
# Save BIC results
write.csv(logl_df, file = paste0(output_file_dir, "summary_LogLikelihood_values.csv"), row.names = FALSE)



#### 8. Check and compare manually extracted BIC/tree BIC results with the automatically extracted values ####
compare_BIC_file = paste0(output_file_dir, "qualityCheck_ML_BIC_comparisons_error.csv")
# Create the compare_BIC_df if it doesn't exist
if (file.exists(compare_BIC_df) == FALSE){
  # Set order for models and datasets
  model_order <- c("PMSF_C20", "PMSF_C60", "PMSF_LG_C20", "PMSF_LG_C60", 
                   "C20", "C60", "LG_C20", "LG_C60", "CF4", "EHO", "EX_EHO",
                   "EX2", "EX3", "GTR20", "JTT", "JTTDCMut", "LG", "LG4M",
                   "mtZOA", "PMB", "Poisson", "rtREV", "UL2", "UL3", "WAG",
                   "ModelFinder")
  dataset_order <- c("Dunn2008", "Philippe2009", "Pick2010", "Philippe2011",
                     "Nosenko2013", "Ryan2013", "Moroz2014", "Borowiec2015",
                     "Chang2015", "Whelan2015", "Whelan2017", "Laumer2018",
                     "Laumer2019")
  # Open the ML tree results and trim the dataframes
  ml_results_file <- grep("maximum_likelihood_results", all_output_files, value = T)
  ml_results_df <- read.table(ml_results_file, header = TRUE, sep = "\t")
  # Fix matrix names to be identical to those in the manual check .xlsx file
  ml_results_df$matrix_name[which(ml_results_df$matrix_name == "Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned")] <- "Dataset10"
  # Trim the ML tree results dataframe
  ml_keep_rows <- sort(intersect(which(ml_results_df$dataset %in% dataset_order), which(ml_results_df$model_code %in% model_order)))
  ml_trimmed_df <- ml_results_df[ml_keep_rows, ]
  # Read in .xlsx file with manual topology check results
  topology_check_x_file <- grep("xls", grep("ML_tree_topology_ManualCheck", all_output_files, value = TRUE), value = TRUE)
  topology_check_x_file <- grep("Summary", topology_check_x_file, value = TRUE, invert = TRUE)
  topology_check_df <- as.data.frame(read_excel(path = topology_check_x_file, sheet = "Topology"))
  # Remove Simion and Hejnol datasets - too computationally intensive to run full ML models
  topology_check_df <- topology_check_df[which(topology_check_df$dataset != "Hejnol2009" & topology_check_df$dataset != "Simion2017"), ]
  # Trim the topology_check_df dataframe
  t_keep_rows <- sort(intersect(which(topology_check_df$dataset %in% dataset_order), which(topology_check_df$model_code %in% model_order)))
  topology_trimmed_df <- topology_check_df[t_keep_rows, ]
  # Update prefix columns
  ml_trimmed_df$prefix <- paste0(ml_trimmed_df$dataset, ".", ml_trimmed_df$matrix_name, ".", ml_trimmed_df$model_code)
  # Trim topology_df to only inclue columns in the ml_trimmed_df
  topology_trimmed_df <- topology_trimmed_df[which(topology_trimmed_df$ML_tree_name %in% ml_trimmed_df$prefix), ]
  # Sort dataframe row order
  topology_trimmed_df <- topology_trimmed_df[order(topology_trimmed_df$dataset, topology_trimmed_df$matrix_name, topology_trimmed_df$model_code), ]
  ml_trimmed_df <- ml_trimmed_df[order(ml_trimmed_df$dataset, ml_trimmed_df$matrix_name, ml_trimmed_df$model_code), ]
  # Combine relevant columns into a single dataframe
  comparison_df <- as.data.frame(cbind(ml_trimmed_df$dataset, ml_trimmed_df$model_code, 
                                       ml_trimmed_df$matrix_name, ml_trimmed_df$prefix,
                                       ml_trimmed_df$best_model_BIC, topology_trimmed_df$model_BIC,
                                       ml_trimmed_df$tree_BIC, topology_trimmed_df$tree_BIC))
  names(comparison_df) <- c("dataset", "model_code",
                            "matrix_name", "prefix",
                            "extracted_model_BIC", "manual_model_BIC",
                            "extracted_tree_BIC", "manual_tree_BIC")
  # Identify which model_BIC and tree_BIC differ in the two files and print results
  differing_model_BIC_rows <- which(signif(as.numeric(comparison_df$extracted_model_BIC, digits = 10)) != signif(as.numeric(comparison_df$manual_model_BIC, digits = 10)) )
  differing_tree_BIC_rows <- which(signif(as.numeric(comparison_df$extracted_tree_BIC, digits = 11)) != signif(as.numeric(comparison_df$manual_tree_BIC, digits = 11)) )
  print(paste0("Number of best model BIC values differing between automatic and manual extraction: ", length(differing_model_BIC_rows)))
  if (length(differing_model_BIC_rows) > 0){print(paste0("Rows differing for model BIC: ", differing_model_BIC_rows))}
  print(paste0("Number of tree BIC values differing between automatic and manual extraction: ", length(differing_tree_BIC_rows)))
  if (length(differing_tree_BIC_rows) > 0){print(paste0("Rows differing for tree BIC: ", differing_tree_BIC_rows))}
  # Output csv if any rows differ
  different_df <- comparison_df[c(differing_model_BIC_rows, differing_tree_BIC_rows), ]
  write.csv(different_df, file = compare_BIC_file, row.names = FALSE)
  # Note: must manually check and correct any rows with differing BIC values (for best model or for trees)
}



#### 9. Collate MAST analyses and constrained tree analyses to calculate BIC ####
# Read in tsv file
mast_df                   <- read.csv(file = grep("04_01_MAST_model_output_SubstitutionModels.csv", all_output_files, value = TRUE), header = TRUE)
# Update mast df
mast_df$tree_topology     <- NA
mast_df2                  <- mast_df[ , c("dataset", "matrix_name", "number_hypothesis_trees", "tree_topology", "model_class", "model_code", "subs_model",
                                          "subs_model_num_params", "mixture_component", "mixture_component_num_params", "rate_num_params",
                                          "state_freq", "state_freq_num_params", "mast_branch_type", "log_likelihood_tree",
                                          "unconstrained_log_likelihood", "num_free_params", "AIC", "AICc", "BIC")]
names(mast_df2)           <-  c("dataset", "matrix_name", "num_trees", "tree_topology", "model_class", "model_code", "model",
                                "subs_model_num_params", "mixture_component", "mixture_component_num_params", "rate_num_params",
                                "state_freq", "state_freq_num_params", "mast_branch_type", "tree_LogL", 
                                "unconstrained_LogL", "num_free_params", "AIC", "AICc", "BIC")
mast_df2$mast_branch_type <- "TR"
# Extract list of constrained tree hypotheses
all_h_trees         <- grep("00_", list.files(hypothesis_tree_dir, recursive = TRUE), value = TRUE, invert = TRUE)
h_tree_iqtree_files <- grep("\\.iqtree", all_h_trees, value = T)
extract_h_trees     <- paste0(hypothesis_tree_dir, grep("ML_H1|ML_H2", h_tree_iqtree_files, value = TRUE))
htree_df            <- as.data.frame(do.call(rbind, lapply(extract_h_trees, extract.hypothesis.tree.parameters)))
htree_df                  <- htree_df[ , c("dataset", "matrix_name", "num_trees", "tree_topology", "model_class", "model_code", 
                                           "subs_model", "subs_model_num_params", "mixture_component", "mixture_component_num_params", "rate_num_params",
                                           "state_freq", "state_freq_num_params", "mast_branch_type", "LogL",
                                           "Unconstrained_LogL", "NumFreeParams", "AIC", "AICc", "BIC")]
names(htree_df)           <-  c("dataset", "matrix_name", "num_trees", "tree_topology", "model_class", "model_code", "model",
                                "subs_model_num_params", "mixture_component", "mixture_component_num_params", "rate_num_params",
                                "state_freq", "state_freq_num_params", "mast_branch_type", "tree_LogL", 
                                "unconstrained_LogL", "num_free_params", "AIC", "AICc", "BIC")
# Concatenate dataframes
check_BIC_df <- rbind(mast_df2, htree_df)
# Order dataframe by dataset and model
check_BIC_df <- check_BIC_df[order(check_BIC_df$dataset, check_BIC_df$matrix_name, check_BIC_df$model_class, 
                                   check_BIC_df$num_trees, check_BIC_df$tree_topology, decreasing = FALSE), ]
rownames(check_BIC_df) <- 1:nrow(check_BIC_df)
# Save dataframe
collated_file <- paste0(output_file_dir, "05_collated_trees_and_MAST.csv")
check_BIC_file <- paste0(output_file_dir, "05_recalculate_BIC.csv")
write.csv(check_BIC_df, file = collated_file)
write.csv(check_BIC_df, file = check_BIC_file)


