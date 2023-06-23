## caitlinch/metazoan-mixtures/code/func_data_analysis.R
# Functions for data analysis and manipulation in R
# Caitlin Cherryh 2023


#### Packages ####
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


extract.number.informative.sites <- function(iq_file){
  ## Extract number of informative sites from IQ-Tree output files
  
  # Get details about the alignment from the filepath name
  dataset = strsplit(basename(iq_file), "\\.")[[1]][1]
  sequence_format = strsplit(basename(iq_file), "\\.")[[1]][3]
  matrix_name <- strsplit(basename(iq_file), "\\.")[[1]][2]
  # Open IQ-Tree file
  iq_lines <- readLines(iq_file)
  # Extract information
  num_sites <- as.numeric(gsub(" ", "", strsplit(strsplit(grep("Input data", iq_lines, value = T), "with")[[1]][2], " ")[[1]][2]))
  num_constant_sites <- as.numeric(gsub(" ", "", strsplit(strsplit(grep("Number of constant sites", iq_lines, value = T), ":")[[1]][2], "\\(")[[1]][1]))
  prop_constant_sites <- num_constant_sites/num_sites
  num_invariant_sites <- as.numeric(gsub(" ", "", strsplit(strsplit(grep("Number of invariant", iq_lines, value = T), ":")[[1]][2], "\\(")[[1]][1]))
  prop_invariant_sites <- num_invariant_sites/num_sites
  num_pis <- as.numeric(gsub(" ", "", strsplit(strsplit(grep("Number of parsimony informative sites:", iq_lines, value = T), ":")[[1]][2], "\\(")[[1]]))
  prop_pis <- num_pis/num_sites
  
  # Collate results into a vector
  op_vector <- c(dataset, matrix_name, sequence_format, num_sites, 
                 num_constant_sites, prop_constant_sites,
                 num_invariant_sites, prop_invariant_sites, 
                 num_pis, prop_pis, iq_file)
  names(op_vector) <- c("dataset", "matrix_name", "sequence_format", "num_sites",
                            "number_constant_sites", "proportion_constant_sites",
                            "number_invariant_sites", "proportion_invariant_sites",
                            "number_informative_sites", "proportion_informative_sites")
  
  # Return the results
  return(op_vector)
}

