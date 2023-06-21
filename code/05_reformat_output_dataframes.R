## caitlinch/metazoan-mixtures/code/05_reformat_output_dataframes.R
# This script reads in csv/tsv files created prior in the pipeline, and reformats them nicely for the manuscript
# Caitlin Cherryh 2023


#### 1. Input parameters ####
## File paths
# output_file_dir         <- Directory for output csvs
# repo_dir                <- Location of caitlinch/metazoan-mixtures github repository

location = "local"
if (location == "local"){
  ## File paths
  output_file_dir         <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/01_output_files/"
  repo_dir                <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"
} else if (location == "dayhoff"){
  ## File paths
  output_file_dir         <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/output/"
  repo_dir                <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/"
  
}



#### 2. Prepare variables, open packages and source functions ####
# Change the default number of digits (so we don't lose the decimal points in the BIC scores)
options(digits = 12)

# Open packages
library(readxl)

# Source files containing functions
source(paste0(repo_dir, "code/func_data_processing.R"))

# List all files in output directory
all_output_files <- paste0(output_file_dir, list.files(output_file_dir))



#### 3. Prepare summary of the manual topology check csv file ####
### Summarise topology results (as percentage of each output topology)
# Read in .xlsx file with manual topology check results
topology_check_x_file <- grep("xls", grep("ML_tree_topology_ManualCheck", all_output_files, value = TRUE), value = TRUE)
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
summary_topology_file <- paste0(output_file_dir, "summary_ML_tree_topology.csv")
write.csv(summary_topology_df, file = summary_topology_file, row.names = FALSE)

### Nicely format output data frames of tree topologies and of sponge topologies
model_order <- c("PMSF_C20", "PMSF_C60", "PMSF_LG_C20", "PMSF_LG_C60", 
                 "C20", "C60", "LG_C20", "LG_C60", "CF4", "EHO", "EX_EHO",
                 "EX2", "EX3", "GTR20", "JTT", "JTTDCMut", "LG", "LG4M",
                 "mtZOA", "PMB", "Poisson", "rtREV", "UL2", "UL3", "WAG",
                 "ModelFinder")
ordered_dataset_ids <- paste0(summary_topology_df$dataset, ".", summary_topology_df$matrix_name)
## For tree topologies:
# Extract tree topologies from dataframe
tree_topology_list <- lapply(ordered_dataset_ids, tree.topology.results, topology_check_df, model_order)
tree_topology_df <- as.data.frame(do.call(cbind, tree_topology_list))
# Format dataframe
names(tree_topology_df) <- ordered_dataset_ids
tree_topology_df$row_names <- c("dataset", "matrix_name", model_order)
tree_topology_df <- tree_topology_df[, c("row_names", ordered_dataset_ids)]
# Output dataframe
tree_topology_file <- paste0(output_file_dir, "all_models_ML_tree_topology.csv")
write.csv(tree_topology_df, file = tree_topology_file, row.names = FALSE)
## For Porifera topologies:
# Extract tree topologies from dataframe
pori_topology_list <- lapply(ordered_dataset_ids, porifera.topology.results, topology_check_df, model_order)
pori_topology_df <- as.data.frame(do.call(cbind, pori_topology_list))
# Format dataframe
names(pori_topology_df) <- ordered_dataset_ids
pori_topology_df$row_names <- c("dataset", "matrix_name", model_order)
pori_topology_df <- pori_topology_df[, c("row_names", ordered_dataset_ids)]
# Output dataframe
pori_topology_file <- paste0(output_file_dir, "all_models_ML_Porifera_topology.csv")
write.csv(pori_topology_df, file = pori_topology_file, row.names = FALSE)



#### 4. Prepare summary of the tree topology test tsv file ####
# Read in tsv file
au_test_file <- grep("tree_topology_test_results.tsv", all_output_files, value = TRUE)
au_test_df <- read.table(file = au_test_file, header = TRUE, sep = "\t")
# Process each dataset one at a time
summary_au_test_list <- lapply(unique(au_test_df$ID), summarise.AU.test.results, au_test_df)
summary_au_test_df <- as.data.frame(do.call(rbind, summary_au_test_list))
# Sort output by year
summary_au_test_df <- summary_au_test_df[order(summary_au_test_df$year, summary_au_test_df$dataset, summary_au_test_df$matrix),]
# Write the output
summary_au_test_file <- paste0(output_file_dir, "summary_au_test_results.csv")
write.csv(summary_au_test_df, file = summary_au_test_file, row.names = FALSE)



#### 5. Prepare summary of the MAST tsv file ####
### Output site ratios from HMM weights
# Read in tsv file
mast_file <- grep("MAST_model_output", all_output_files, value = TRUE)
mast_df <- read.table(file = mast_file, header = TRUE, sep = "\t")
# Separate into HMMster and phyloHMM dfs
phylohmm_results_df <- mast_df[mast_df$analysis_type == "phyloHMM",]
hmmster_results_df <- mast_df[mast_df$analysis_type == "HMMster",]
# Process each dataset one at a time
phylohmm_summary_df <- as.data.frame(do.call(rbind, lapply(1:nrow(phylohmm_results_df), summarise.HMM.results, hmm_df = phylohmm_results_df)))
hmmster_summary_df <- as.data.frame(do.call(rbind, lapply(1:nrow(hmmster_results_df), summarise.HMM.results, hmm_df = hmmster_results_df)))
# Add branch length options
phylohmm_summary_df$branch_length_option <- "TR"
hmmster_summary_df$branch_length_option <- "T"
# Sort output by year
phylohmm_summary_df <- phylohmm_summary_df[order(phylohmm_summary_df$year, phylohmm_summary_df$dataset, phylohmm_summary_df$matrix_name), ]
hmmster_summary_df <- hmmster_summary_df[order(hmmster_summary_df$year, hmmster_summary_df$dataset, hmmster_summary_df$matrix_name), ]
# Relabel row numbers
rownames(phylohmm_summary_df) <- 1:nrow(phylohmm_summary_df)
rownames(hmmster_summary_df) <- 1:nrow(hmmster_summary_df)
# Write the output for phyloHMM
phylohmm_summary_test_file <- paste0(output_file_dir, "summary_phyloHMM_HMM_results.csv")
write.csv(phylohmm_summary_df, file = phylohmm_summary_test_file, row.names = FALSE)
# Write the output for HMMster
hmmster_summary_test_file <- paste0(output_file_dir, "summary_HMMster_HMM_results.csv")
write.csv(hmmster_summary_df, file = hmmster_summary_test_file, row.names = FALSE)

### Output tree weights
# Process each dataset one at a time
phylohmm_summary_df <- as.data.frame(do.call(rbind, lapply(1:nrow(phylohmm_results_df), summarise.tree.weights, tw_df = phylohmm_results_df)))
hmmster_summary_df <- as.data.frame(do.call(rbind, lapply(1:nrow(hmmster_results_df), summarise.tree.weights, tw_df = hmmster_results_df)))
# Add branch length options
phylohmm_summary_df$branch_length_option <- "TR"
hmmster_summary_df$branch_length_option <- "T"
# Sort output by year
phylohmm_summary_df <- phylohmm_summary_df[order(phylohmm_summary_df$year, phylohmm_summary_df$dataset, phylohmm_summary_df$matrix_name), ]
hmmster_summary_df <- hmmster_summary_df[order(hmmster_summary_df$year, hmmster_summary_df$dataset, hmmster_summary_df$matrix_name), ]
# Relabel row numbers
rownames(phylohmm_summary_df) <- 1:nrow(phylohmm_summary_df)
rownames(hmmster_summary_df) <- 1:nrow(hmmster_summary_df)
# Write the output for phyloHMM
phylohmm_summary_test_file <- paste0(output_file_dir, "summary_phyloHMM_treeWeight_results.csv")
write.csv(phylohmm_summary_df, file = phylohmm_summary_test_file, row.names = FALSE)
# Write the output for HMMster
hmmster_summary_test_file <- paste0(output_file_dir, "summary_HMMster_treeWeight_results.csv")
write.csv(hmmster_summary_df, file = hmmster_summary_test_file, row.names = FALSE)



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


