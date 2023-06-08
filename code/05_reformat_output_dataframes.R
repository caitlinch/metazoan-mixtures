# metazoan-mixtures/code/05_reformat_output_dataframes.R
## This script reads in tsv files created prior in the pipeline, and reformats them nicely for the manuscript
# Caitlin Cherryh, 2023



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
# Source files containing functions
source(paste0(repo_dir, "code/func_data_processing.R"))

# List all files in output directory
all_output_files <- paste0(output_file_dir, list.files(output_file_dir))



#### 3. Prepare summary of the manual topology check csv file ####
# Read in csv file
topology_check_file <- grep("ML_tree_topology_ManualCheck", all_output_files, value = TRUE)
topology_check_df <- read.csv(topology_check_file, stringsAsFactors = FALSE)



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




