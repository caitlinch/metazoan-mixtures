## caitlinch/metazoan-mixtures/func_pmsf_trees.R
# Caitlin Cherryh 2022

library(ape) # Functions: read.nexus.data
library(phylotools) # Functions: read.fasta, dat2phylip

# Functions for estimating maximum likelihood trees and constrained maximum likelihood trees with PMSF and CAT-PMSF models



#### Functions for estimating trees with PMSF model in IQ-Tree ####
estimate.PMSF.tree.wrapper <- function(row_id, pmsf_parameter_dataframe, run.iqtree = FALSE){
  # Function to apply the estimate.PMSF.tree function to multiple datasets
  
  # Take one row from the dataframe
  row <- pmsf_parameter_dataframe[row_id, ]
  
  # Apply the PMSF function
  row_pmsf_output <- estimate.PMSF.tree(alignment_path = row$alignment_file, alignment_prefix = row$prefix, simple_model = row$guide_tree_model,
                                        iqtree_path = row$iqtree_path, num_threads = row$iqtree_num_threads, num_bootstraps = row$iqtree_num_bootstraps,
                                        pmsf_dir = row$pmsf_dir, run.iqtree = FALSE)
  
  # Return the output
  return(row_pmsf_output)
}


estimate.PMSF.tree <- function(alignment_path, alignment_prefix, simple_model, iqtree_path, num_threads, num_bootstraps, pmsf_dir, run.iqtree = FALSE){
  # Function to estimate a ML tree using the PMSF (posterior mean site frequency) model, from start to finish
  
  # Change working location to output directory
  setwd(pmsf_dir)
  
  # 1. Estimate guide tree under simple model
  guide_tree_command <- estimate.guide.tree(alignment_path, alignment_prefix, simple_model, iqtree_path, num_threads, run.iqtree = FALSE)[[1]]
  guide_tree_prefix <- estimate.guide.tree(alignment_path, alignment_prefix, simple_model, iqtree_path, num_threads, run.iqtree = FALSE)[[2]]
  # List all files in the PMSF directory
  all_files_pmsf_dir <- list.files(pmsf_dir)
  # Find the guide tree file
  guide_tree_check <- grep("treefile", grep(guide_tree_prefix, all_files_pmsf_dir, value = T), value = T)
  if (length(guide_tree_check) == 0){
    # No guide tree file: return NA
    guide_tree_path = NA
  } else if (length(guide_tree_check) > 0){
    # Guide tree file exists: return guide tree file
    guide_tree_path <- paste0(pmsf_dir, guide_tree_check)
  }
  
  # 2. Perform the first phase of the PMSF model: estimate mixture model parameters given the guide tree and infer site-specific 
  #   frequency profile (printed to .sitefreq file)
  sitefreq_command <- output.site.frequency.file(alignment_path, guide_tree_path, alignment_prefix, simple_model, iqtree_path, num_threads, run.iqtree = FALSE)[[1]]
  sitefreq_prefix <- output.site.frequency.file(alignment_path, guide_tree_path, alignment_prefix, simple_model, iqtree_path, num_threads, run.iqtree = FALSE)[[2]]
  # List all files in the PMSF directory
  all_files_pmsf_dir <- list.files(pmsf_dir)
  # Find the sitefreq file
  sitefreq_check <- grep("treefile", grep(sitefreq_prefix, all_files_pmsf_dir, value = T), value = T)
  if (length(sitefreq_check) == 0){
    # No guide tree file: return NA
    sitefreq_path = NA
  } else if (length(sitefreq_check) == 1){
    # Guide tree file exists: return guide tree file
    sitefreq_path <- paste0(pmsf_dir, sitefreq_check)
  }
  
  # 3. Perform the second phase of the PMSF model: conduct typical analysis using the inferred frequency model (instead of the mixture model) 
  #   to save RAM and running time.
  pmsf_command <- estimate.tree.with.inferred.PMSF.model(alignment_path, sitefreq_path, alignment_prefix, simple_model, iqtree_path, num_threads, num_bootstraps, run.iqtree = FALSE)[[1]]
  pmsf_prefix <- estimate.tree.with.inferred.PMSF.model(alignment_path, sitefreq_path, alignment_prefix, simple_model, iqtree_path, num_threads, num_bootstraps, run.iqtree = FALSE)[[2]]
  # List all files in the PMSF directory
  all_files_pmsf_dir <- list.files(pmsf_dir)
  # Find the pmsf tree files
  # Find the sitefreq file
  pmsf_check <- grep("iqtree", grep(pmsf_prefix, all_files_pmsf_dir, value = T), value = T)
  if (length(sitefreq_check) == 0){
    # No PMSF tree file: return NA for output files with PMSF prefix
    pmsf_treefile <- NA
    pmsf_iqfile <- NA
    pmsf_logfile <- NA
  } else if (length(sitefreq_check) == 1){
    # PMSF tree file exists: return output files with PMSF prefix
    pmsf_treefile <- paste0(pmsf_dir, grep("treefile", grep(pmsf_prefix, all_files_pmsf_dir, value = T), value = T))
    pmsf_iqfile <- paste0(pmsf_dir, grep("iqtree", grep(pmsf_prefix, all_files_pmsf_dir, value = T), value = T))
    pmsf_logfile <- paste0(pmsf_dir, grep("log", grep(pmsf_prefix, all_files_pmsf_dir, value = T), value = T))
  }
  
  # Create a little output vector of all the information
  output_vector <- c(alignment_prefix, alignment_path, "PMSF", simple_model, num_threads, num_bootstraps, guide_tree_command, guide_tree_path, 
                     sitefreq_command, sitefreq_path, pmsf_command, pmsf_iqfile, pmsf_logfile, pmsf_treefile, pmsf_dir)
  names(output_vector) =c("prefix", "alignment_path", "model_code", "input_model", "num_threads", "num_bootstraps", "IQTree_command_1", "guide_tree_path", 
                          "IQTree_command_2", "site_frequencies_path", "IQTree_command_3", "pmsf_iqtree_file", "pmsf_log_file", "pmsf_tree_file", "pmsf_directory")
  # Return the output vector
  return(output_vector)
}


find.pmsf.files <- function(pmsf_prefix, pmsf_dir){
  # Quick function to return the PMSF output files for a single run, using the alignment prefix
  
  # List all files in the PMSF directory
  all_files_pmsf_dir <- list.files(pmsf_dir)
  # Find the pmsf tree files
  # Find the sitefreq file
  pmsf_check <- grep("iqtree", grep(pmsf_prefix, all_files_pmsf_dir, value = T), value = T)
  if (length(sitefreq_check) == 0){
    # No PMSF tree file: return NA for output files with PMSF prefix
    pmsf_treefile <- NA
    pmsf_iqfile <- NA
    pmsf_logfile <- NA
  } else if (length(sitefreq_check) == 1){
    # PMSF tree file exists: return output files with PMSF prefix
    pmsf_treefile <- paste0(pmsf_dir, grep("treefile", grep(pmsf_prefix, all_files_pmsf_dir, value = T), value = T))
    pmsf_iqfile <- paste0(pmsf_dir, grep("iqtree", grep(pmsf_prefix, all_files_pmsf_dir, value = T), value = T))
    pmsf_logfile <- paste0(pmsf_dir, grep("log", grep(pmsf_prefix, all_files_pmsf_dir, value = T), value = T))
  }
  
  # Return the three output files
  return(c(pmsf_treefile, pmsf_iqfile, pmsf_logfile))
}


estimate.guide.tree.wrapper <- function(row_id, pmsf_parameter_dataframe, run.iqtree = FALSE){
  # Function to apply the estimate.guide.tree function to multiple datasets
  
  # Take one row from the dataframe
  row <- pmsf_parameter_dataframe[row_id, ]
  
  # Run the estimate.guide.tree function
  output_vector <- estimate.guide.tree(alignment_path = row$alignment_file, alignment_prefix = row$prefix, simple_model = row$guide_tree_model, 
                                       iqtree_path = row$iqtree_path, num_threads = row$iqtree_num_threads, run.iqtree = FALSE)
  
  # Return the output
  return(output_vector)
}


estimate.guide.tree <- function(alignment_path, alignment_prefix, simple_model, iqtree_path, num_threads, run.iqtree = FALSE){
  # Function to estimate a guide tree for the PMSF model
  # IQ-Tree command: iqtree -s <alignment> -m 'LG+C20+F+G' -pre guidetree
  
  # Create the prefix for a guide tree
  guide_prefix <- paste0(alignment_prefix, ".guidetree")
  
  # Assemble the command to estimate a guide tree
  iqtree_command <- paste0(iqtree_path, " -s ", alignment_path, " -m ", simple_model, " -nt ", num_threads, " -pre ", guide_prefix)
  
  # Run IQ-Tree2 (if "run.iqtree" == TRUE)
  if (run.iqtree == TRUE){
    system(iqtree_command)
  }
  
  # Return the IQ-Tree command to estimate a guide tree
  return(c(iqtree_command, guide_prefix))
}


output.site.frequency.file.wrapper <- function(row_id, pmsf_parameter_dataframe, run.iqtree = FALSE){
  # Function to apply the output.site.frequency.file function to multiple datasets
  
  # Take one row from the dataframe
  row <- pmsf_parameter_dataframe[row_id, ]
  
  # Run the estimate.guide.tree function
  output_vector <- output.site.frequency.file(alignment_path = row$alignment_file, guide_tree_path = row$guide_tree_path, alignment_prefix = row$prefix, 
                                              simple_model = row$guide_tree_model, iqtree_path = row$iqtree_path, num_threads = row$iqtree_num_threads, run.iqtree = FALSE)
  
  # Return the output
  return(output_vector)
}


output.site.frequency.file <- function(alignment_path, guide_tree_path, alignment_prefix, simple_model, iqtree_path, num_threads, run.iqtree = FALSE){
  # Function to estimate a site frequency file for the PMSF model, given a guide tree
  # IQ-Tree command: iqtree -s <alignment> -m 'LG+C20+F+G' -ft <guide_tree> -n 0 -pre ssfp
  
  # Create the prefix for a sitefreq file (site-specific frequency profile or ssfp)
  ssfp_prefix <- paste0(alignment_prefix, ".ssfp")
  
  # Assemble the command to estimate a guide tree
  iqtree_command <- paste0(iqtree_path, " -s ", alignment_path, " -m ", simple_model, " -ft ", guide_tree_path, " -n 0 -nt ", num_threads, " -pre ", ssfp_prefix)
  
  # Run IQ-Tree2 (if "run.iqtree" == TRUE)
  if (run.iqtree == TRUE){
    system(iqtree_command)
  }
  
  # Return the IQ-Tree command to estimate a guide tree
  return(c(iqtree_command, ssfp_prefix))
}


estimate.tree.with.inferred.PMSF.model.wrapper <- function(row_id, pmsf_parameter_dataframe, run.iqtree = FALSE){
  # Function to apply the estimate.tree.with.inferred.PMSF.model function to multiple datasets
  
  # Take one row from the dataframe
  row <- pmsf_parameter_dataframe[row_id, ]
  
  # Run the estimate.guide.tree function
  output_vector <- estimate.tree.with.inferred.PMSF.model(alignment_path = row$alignment_file, sitefreq_path = row$site_frequencies_path, alignment_prefix = row$prefix, simple_model = row$guide_tree_model, 
                                                          iqtree_path = row$iqtree_path, num_threads = row$iqtree_num_threads, num_bootstraps = row$iqtree_num_bootstraps, run.iqtree = FALSE)
  
  # Return the output
  return(output_vector)
}


estimate.tree.with.inferred.PMSF.model <- function(alignment_path, sitefreq_path, alignment_prefix, simple_model, iqtree_path, num_threads, num_bootstraps, run.iqtree = FALSE){
  # Function to estimate a tree with the PMSF model, using the site frequency file inferred from a guide tree
  # IQ-Tree command: iqtree -s <alignment> -m LG+C20+F+G -fs <file.sitefreq> -b 100
  
  # Create the prefix for the final output tree 
  pmsf_prefix <- paste0(alignment_prefix, ".complete")
  
  # Assemble the command to estimate a guide tree
  iqtree_command <- paste0(iqtree_path, " -s ", alignment_path, " -m ", simple_model, " -fs ", sitefreq_path, " -b ", num_bootstraps, " -nt ", num_threads," -pre ", pmsf_prefix)
  
  # Run IQ-Tree2 (if "run.iqtree" == TRUE)
  if (run.iqtree == TRUE){
    system(iqtree_command)
  }
  
  # Return the IQ-Tree command to estimate a guide tree
  return(c(iqtree_command, pmsf_prefix))
}


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

