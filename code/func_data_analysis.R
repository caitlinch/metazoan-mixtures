## caitlinch/metazoan-mixtures/func_data_analysis.R
# Caitlin Cherryh 2022

# Functions for data analysis and manipulation in R

library(ape)
library(phylotools)

#### Functions to extract information about alignments ####
matrix.dimensions <- function(alignment_path){
  # Small function to  check the number of dimensions in an alignment
  
  # Print the file path
  print(alignment_path)
  
  # Get details about the alignment from the filepath name
  dataset = strsplit(basename(alignment_path), "\\.")[[1]][1]
  sequence_format = strsplit(basename(alignment_path), "\\.")[[1]][3]
  matrix_name <- strsplit(basename(alignment_path), "\\.")[[1]][2]
  
  # Get the file type
  suffix <- tail(strsplit(alignment_path, "\\.")[[1]], 1)
  
  if (suffix == "fa" | suffix == "fas" | suffix == "fasta"){
    # For fasta files
    f <- read.fasta(alignment_path)
    num_taxa <- length(f$seq.name)
    num_sites <- nchar(f$seq.text[[1]])
  } else if (suffix == "nex" | suffix == "nexus"){
    # For nexus files
    n <- read.nexus.data(alignment_path)
    num_taxa <- length(names(n))
    num_sites <- length(n[[1]])
  } else if (suffix == "phy" | suffix == "phylip"){
    # For phylip files
    p <- read.phylip(alignment_path)
    num_taxa <- length(p$seq.name)
    num_sites <- nchar(p$seq.text[[1]])
  }
  
  # Collate results into a vector
  op_vector <- c(dataset, matrix_name, sequence_format, num_taxa, num_sites, alignment_path)
  
  # Return the results
  return(op_vector)
}

