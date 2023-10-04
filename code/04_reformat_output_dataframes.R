## caitlinch/metazoan-mixtures/code/05_reformat_output_dataframes.R
# This script reads in csv/tsv files created prior in the pipeline, and reformats them nicely for the manuscript
# Caitlin Cherryh 2023


#### 1. Input parameters ####
## File paths
# output_file_dir         <- Directory for output csvs
# repo_dir                <- Location of caitlinch/metazoan-mixtures github repository

## File paths
output_file_dir         <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/01_output_files/"
repo_dir                <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"



#### 2. Prepare variables, open packages and source functions ####
# Change the default number of digits (so we don't lose the decimal points in the BIC scores)
options(digits = 12)

# Open packages
library(readxl)

# Source files containing functions
source(paste0(repo_dir, "code/func_data_processing.R"))

# List all files in output directory
all_output_files <- paste0(output_file_dir, list.files(output_file_dir))
# Remove any files with all 5 trees - only want to look at the output for the 2 tree model
all_output_files <- grep("5trees", all_output_files, value = TRUE, invert = TRUE)



#### 3. Prepare summary of the manual topology check csv file ####
### Summarise topology results (as percentage of each output topology)
summary_topology_file <- paste0(output_file_dir, "summary_ML_tree_topology.csv")
if (file.exists(summary_topology_file) == FALSE){
  # Read in .xlsx file with manual topology check results
  topology_check_x_file <- grep("xls", grep("ML_tree_topology_ManualCheck", all_output_files, value = TRUE), value = TRUE)
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
  # Summarise results for each dataset as a percentage
  dataset_ids <- unique(paste0(topology_check_df$dataset, ".", topology_check_df$matrix_name))
  summary_topology_list <- lapply(dataset_ids, summarise.topology.results, topology_check_df, 
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
                 "EX2", "EX3", "GTR20", "JTT", "JTTDCMut", "LG", "LG4M",
                 "mtZOA", "PMB", "Poisson", "rtREV", "UL2", "UL3", "WAG",
                 "ModelFinder")

# Order dataframe by dataset/matrix
ordered_dataset_ids <- paste0(summary_topology_df$dataset, ".", summary_topology_df$matrix_name)
## For tree topologies:
tree_topology_file <- paste0(output_file_dir, "all_models_ML_tree_topology.csv")
if (file.exists(tree_topology_file) == FALSE){
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
# Process each dataset one at a time
summary_au_test_df <- as.data.frame(do.call(rbind, lapply(unique(tree_test_df$ID), summarise.AU.test.results, tree_test_df)))
# Add new column for model class
summary_au_test_df$model_class <- factor(summary_au_test_df$best_model_code,
                                         levels = c("LG_C60", "C60", "PMSF_C60", "PMSF_LG_C60", "LG4M", "UL3"),
                                         labels = c("CXX", "CXX", "PMSF", "PMSF", "Other", "Other"),
                                         ordered = TRUE)
# Sort output by year
summary_au_test_df <- summary_au_test_df[order(summary_au_test_df$year, summary_au_test_df$dataset, summary_au_test_df$matrix),]
# Reorder columns
summary_au_test_df <- summary_au_test_df[, c("dataset", "matrix", "model_class", "best_model_code", "topology_test", "tree_1", "tree_2", "tree_3", "tree_4", "tree_5", "year")]
# Remove columns consisting only of NA
summary_au_test_df <- Filter(function(x)!all(is.na(x)), summary_au_test_df)
# Write the output
summary_au_test_file <- paste0(output_file_dir, "summary_au_test_results.csv")
write.csv(summary_au_test_df, file = summary_au_test_file, row.names = FALSE)



#### 5. Prepare summary of the expected likelihood weights using the tree topology test tsv file ####
# Read in csv file
tree_test_df <- read.csv(file = grep("tree_topology_test_results.csv", all_output_files, value = TRUE), header = TRUE)
# Process each dataset one at a time
summary_elw_test_df <- as.data.frame(do.call(rbind, lapply(unique(tree_test_df$ID), summarise.eLW, tree_test_df)))
# Add new column for model class
summary_elw_test_df$model_class <- factor(summary_elw_test_df$best_model_code,
                                          levels = c("LG_C60", "C60", "PMSF_C60", "PMSF_LG_C60", "LG4M", "UL3"),
                                          labels = c("CXX", "CXX", "PMSF", "PMSF", "Other", "Other"),
                                          ordered = TRUE)
# Sort output by year
summary_elw_test_df <- summary_elw_test_df[order(summary_elw_test_df$year, summary_elw_test_df$dataset, summary_elw_test_df$matrix),]
# Reorder columns
summary_elw_test_df <- summary_elw_test_df[, c("dataset", "matrix", "model_class", "best_model_code", "topology_test", "tree_1", "tree_2", "tree_3", "tree_4", "tree_5", "year")]
# Remove columns consisting only of NA
summary_elw_test_df <- Filter(function(x)!all(is.na(x)), summary_elw_test_df)
# Write the output
summary_elw_test_file <- paste0(output_file_dir, "summary_elw_results.csv")
write.csv(summary_elw_test_df, file = summary_elw_test_file, row.names = FALSE)



#### 5. Prepare summary of the MAST tsv file ####
### Output site ratios from HMM weights
# Read in tsv file
mast_df <- read.csv(file = grep("MAST_model_output", all_output_files, value = TRUE), header = TRUE)
# Output tree weights
# Process each dataset one at a time
mast_summary_df <- as.data.frame(do.call(rbind, lapply(1:nrow(mast_df), summarise.tree.weights, tw_df = mast_df)))
# Sort output by year
mast_summary_df <- mast_summary_df[order(mast_summary_df$year, mast_summary_df$dataset, mast_summary_df$matrix_name), ]
# Reorganise columns
mast_summary_df <- mast_summary_df[,c("dataset", "matrix_name", "model_class", "model_code", 
                                      "tree_1_tree_weight", "tree_2_tree_weight", "tree_3_tree_weight", 
                                      "tree_4_tree_weight", "tree_5_tree_weight",
                                      "mast_branch_type", "minimum_branch_length", "number_hypothesis_trees", "year")]
# Remove columns consisting only of NA
mast_summary_df <- Filter(function(x)!all(is.na(x)), mast_summary_df)
# Write the output for MAST tree weights
write.csv(mast_summary_df, file = paste0(output_file_dir, "summary_MAST_treeWeight_results.csv"), row.names = FALSE)



#### 6. Check and compare manually extracted BIC/tree BIC results with the automatically extracted values ####
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
differing_model_BIC_rows <- which(comparison_df$extracted_model_BIC != comparison_df$manual_model_BIC)
differing_tree_BIC_rows <- which(comparison_df$extracted_tree_BIC != comparison_df$manual_tree_BIC)
print(paste0("Number of best model BIC values differing between automatic and manual extraction: ", length(differing_model_BIC_rows)))
if (length(differing_model_BIC_rows) > 0){print(paste0("Rows differing for model BIC: ", differing_model_BIC_rows))}
print(paste0("Number of tree BIC values differing between automatic and manual extraction: ", length(differing_tree_BIC_rows)))
if (length(differing_tree_BIC_rows) > 0){print(paste0("Rows differing for tree BIC: ", differing_tree_BIC_rows))}
# Manually check and correct any rows with differing BIC values (for best model or for trees)


