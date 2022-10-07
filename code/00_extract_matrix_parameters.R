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
library(phylotools)
library(TreeTools)
library(phangorn)

# Source functions
source(paste0(repo_dir, "code/func_data_analysis.R"))



#### 3. Prepare alignment paths ####
# Extract the list of all files from the folder containing alignments/models/partition files
all_files <- list.files(alignment_dir)
if (length(all_files) > 0){
  all_files <- paste0(alignment_dir, all_files)
}
# Extract the list of alignments (i.e. files that contain the word "alignment")
all_alignments <- grep("\\.alignment\\.", all_files, value = T)

# Two alignment paths need to be corrected - Dunn2008 and Hejnol2009
# To correct Dunn 2008 (issue in reading file - read in as phyDat and write out as fasta file):
corrected_dunn_file <- grep("FixedNames", grep("Dunn2008", all_alignments, value = T), value = T)
if (file.exists(corrected_dunn_file) == FALSE){
  # Get name of original alignment file
  dunn_al <- grep("Original", grep("Dunn2008", all_alignments, value = T), value = T)
  # Open alignment as phyDat
  dunn_data <- ReadAsPhyDat(dunn_al)
  # Write alignment as fasta file
  write.phyDat(dunn_data, file = corrected_dunn_file, format = "fasta", colsep = "") 
}
# To correct Dunn 2008 (issue in reading file - read in as phyDat and write out as fasta file):
corrected_hejnol_file <- grep("Original", grep("Hejnol2009", all_alignments, value = T), value = T)
if (file.exists(corrected_hejnol_file) == FALSE){
  # Get name of original alignment file
  hejnol_al <- grep("FixedNames", grep("Hejnol2009", all_alignments, value = T), value = T)
  # Open alignment as phyDat
  hejnol_data <- ReadAsPhyDat(hejnol_al)
  # Write alignment as fasta file
  write.phyDat(hejnol_data, file = corrected_hejnol_file, format = "fasta", colsep = "")
}

# Re-extract the list of alignments and remove the files containing the word "Original"
# Extract the list of all files from the folder containing alignments/models/partition files
all_files <- list.files(alignment_dir)
if (length(all_files) > 0){
  all_files <- paste0(alignment_dir, all_files)
}
# Extract the list of alignments (i.e. files that contain the word "alignment")
all_alignments <- grep("\\.alignment\\.", all_files, value = T)
all_alignments <- grep("Original", all_alignments, value = T, invert = T)



#### 4. Extract matrix dimensions ####
# Extract information about each alignment
dimension_list <- lapply(all_alignments, matrix.dimensions)
dimension_df <- as.data.frame(do.call(rbind, dimension_list))
names(dimension_df) <- c("dataset", "matrix_name", "sequence_format", "num_taxa", "num_sites", "alignment_path")

# Save the dataframe
df_path <- paste0(output_dir, "alignment_dimensions.csv")
write.csv(dimension_df, file = df_path, row.names = FALSE)

