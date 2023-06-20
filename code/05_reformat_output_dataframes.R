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



#### 5. Prepare summary of the phyloHMM tsv file ####
# Read in tsv file
phyloHMM_file <- grep("phyloHMM", all_output_files, value = TRUE)
phyloHMM_df <- read.table(file = phyloHMM_file, header = TRUE, sep = "\t")
# Process each dataset one at a time

# Sort output by year

# Write the output



#### 6. Check and compare manually extracted BIC/tree BIC results with the automatically extracted values ####





