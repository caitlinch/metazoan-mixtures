## caitlinch/metazoan-mixtures/func_pmsf_trees.R
# Caitlin Cherryh 2022

library(ape) # Functions: read.nexus.data
library(phylotools) # Functions: read.fasta, dat2phylip

# Functions for estimating maximum likelihood trees and constrained maximum likelihood trees with PMSF and CAT-PMSF models



#### Functions for estimating trees with PMSF model in IQ-Tree ####



#### Functions for estimating trees with the CAT-PMSF model ####



### Functions for managing and manipulating alignments ####
convert.to.phylip <- function(alignment_path, sequence_format = "AA"){
  ## Convert nexus or fasta files to phylip
  
  # Check what the file format of the alignment is
  suffix <- tail(strsplit(basename(alignment_path), "\\.")[[1]],1)
  
  if (suffix == "phy" | suffix == "phylip"){
    ## Alignment is already a phylip file
    # Return the path to the existing phylip file
    phylip_alignment_file_path <- alignment_path
  } else if (suffix == "fasta"| suffix == "fas" | suffix == "fa" | suffix == "fna" | suffix == "ffn" | suffix == "faa" | suffix == "frn"){
    ## Alignment is a fasta file format - convert to phylip
    # Create new file name for output phylip file
    phylip_alignment_file_path <- paste0(alignment_path, ".phylip")
    # Read in the fasta file
    f_data <- read.fasta(file = alignment_path)
    # Write out the fasta data as a phylip file
    dat2phylip(dat = f_data, outfile = phylip_alignment_file_path)
  } else if (suffix == "nex" | suffix == "nexus" | suffix == "nxs"){
    ## Alignment is a nexus file - convert to phylip
    # Create new file name for output files
    phylip_alignment_file_path <- paste0(alignment_path, ".phylip")
    fasta_alignment_file_path <- paste0(alignment_path, ".fasta")
    # Read in the nexus file
    n_data <- read.nexus.data(alignment_path)
    # Change nexus data to AA sequences
    n_AA_seqs <- as.AAbin(n_data)
    # Write data as a fasta file
    write.FASTA(n_AA_seqs, file = fasta_alignment_file_path)
    # Translate the fasta file to a phylip file
    f_data <- read.fasta(file = fasta_alignment_file_path)
    dat2phylip(dat = f_data, outfile = phylip_alignment_file_path)
  }
  
  # Return the path to the phylip alignment file
  return(phylip_alignment_file_path)
}

