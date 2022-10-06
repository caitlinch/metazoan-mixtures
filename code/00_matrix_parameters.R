# metazoan-mixtures/code/00_matrix_parameters.R
## This script extracts the number of taxa and number of sites from a list of alignments
# Caitlin Cherryh, 2022

#### 1. Input parameters ####
## Specify parameters:
# alignment_dir       <- Directory containing alignments for all data sets
#                        Alignments have the naming convention dataset.matrix_name.sequence_type.fa
#                        E.g. Cherryh2022.all_taxa.aa.fa
# output_dir          <- Directory for IQ-Tree output (trees and tree mixtures)
# repo_dir            <- Location of caitlinch/metazoan-mixtures github repository

location = "local"
if (location == "local"){
  alignment_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/01_Data_all/"
  output_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/"
  repo_dir <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"
} else if (location == "soma"){
  alignment_dir <- "/data/caitlin/metazoan-mixtures/data_all/"
  output_dir <- "/data/caitlin/metazoan-mixtures/output/"
  repo_dir <- "/data/caitlin/metazoan-mixtures/"
} 



#### 2. Open libraries ####
# Open packages
library(ape)
library(phylotools)

# Source functions
source(paste0(repo_dir, "code/func_data_analysis.R"))



#### 3. Extract matrix dimensions ####
# Extract the list of all files from the folder containing alignments/models/partition files
all_files <- list.files(alignment_dir)
if (length(all_files) > 0){
  all_files <- paste0(alignment_dir, all_files)
}
# Extract the list of alignments (i.e. files that contain the word "alignment")
all_alignments <- grep("\\.alignment\\.", all_files, value = T)

# Extract information about each alignment
dimension_list <- lapply(all_alignments, matrix.dimensions)
dimension_df <- do.call(rbind, dimension_list)
names(dimension_df) <- c("dataset", "matrix_name", "sequence_format", "num_taxa", "num_sites", "alignment_path")

# Save the dataframe
df_path <- paste0(output_dir, "alignment_dimensions.csv")
write.csv(dimension_df, file = df_path, row.names = FALSE)
