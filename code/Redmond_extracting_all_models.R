# File to extract all models used in the Redmond and McLysaght paper (from the partition files!)
# Paper: https://www.nature.com/articles/s41467-021-22074-7
# Data repository: https://doi.org/10.6084/m9.figshare.12746972

# Input parameters
# supp_dir <- directory containing supplementary data from Redmond paper(downloaded from repository)

supp_dir <- "/Users/caitlin/Downloads/RedmondMcLysaght2021NatCommsData" 
mandir <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"



#### 1. Prepare variables, packages, and functions
source(paste0(main_dir, "code/func_data_processing.R"))



#### 2. Extract vector of all models used in Redmond 2021 ####
# List all files in the supplementary data
all_files <- list.files(supp_dir, recursive = T, full.names = T)
# Extract files for the three sponge datasets (REA, WEA15 and WEA17)
sponge_files <- grep("REA|WEA15|WEA17", all_files, value = T)
# Extract partition files for sponge datasets
sponge_partition_files <- grep("partition", sponge_files, value = T)
# Remove recoded analyses
sponge_partition_files <- grep("RL1|RL2", sponge_partition_files, value = T, invert = T)
# Extract models from all partition files
sponge_partition_model_list <- lapply(sponge_partition_files, extract.partition.models)
# Combine list of partition model dataframes into a single dataframe
models <- sort(unique(unlist(sponge_partition_model_list)))

