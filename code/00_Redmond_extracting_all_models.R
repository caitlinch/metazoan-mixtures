## caitlinch/metazoan-mixtures/code/00_Redmond_extracting_all_models.R
# This script extracts the list of models used for tree estimation in a precious study, 
# Caitlin Cherryh 2023

# File to extract all models used in the Redmond and McLysaght paper (from the partition files!)
# Paper: https://www.nature.com/articles/s41467-021-22074-7
# Data repository: https://doi.org/10.6084/m9.figshare.12746972

#### 1. Input parameters ####
# supp_dir <- Directory containing supplementary data from Redmond and McLysaght (2021) paper 
#             Downloaded from repository here: https://doi.org/10.6084/m9.figshare.12746972
# main_dir <- Directory for GitHub repository for this project (caitlinch/metazoan-mixtures)

supp_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/RedmondMcLysaght2021NatCommsData" 
main_dir <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"



#### 2. Prepare variables, packages, and functions
source(paste0(main_dir, "code/func_data_processing.R"))



#### 3. Extract vector of all models used in Redmond and McLysaght 2021 ####
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

# Ensure the output folder exists (and create it if it doesn't)
output_dir <- paste0(main_dir, "data/")
if (dir.exists(output_dir) == FALSE){
  dir.create(output_dir)
}

# Save the models as a text file
output_path <- paste0(output_dir, "Redmond.McLysaght2021_all_models.txt")
write(models, file = output_path)


