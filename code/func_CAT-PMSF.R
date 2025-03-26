## caitlinch/metazoan-mixtures/code/func_CAT-PMSF.R
# Functions for estimating CAT-PMSF models using PhyloBayes, and inferring trees with those models
# Caitlin Cherryh 2025

library(ape)

infer.phylobayes.profile.wrapper <- function(i, alignment_df, params, profile_topologies){
  # Generate command lines to infer compositional profiles in PhyloBayes for a
  # given dataset

  alignment_details <- alignment_df[i, ]

  # Need to generate 2 profiles (under CTEN and PORI fixed topologies) with 2 chains each:
  # TODO Add profile topologies to four_profiles
  dataset_details <- c(
    "dataset" = alignment_details$dataset,
    "matrix_name" = alignment_details$matrix_name,
    "alignment" = alignment_details$full_path,
    "threads" = params$phylobayes_threads,
    "pb_mpi_dir" = params$phylobayes
  )
  PB_profiles <- list(
    CTEN_profile = c(
      "fixed_tree" = "CTEN",
      "tree" = grep("constraint_tree_1",
                    grep(alignment_details$matrix_name,
                         grep(alignment_details$dataset, profile_topologies, value = T),
                         value = T),
                    value = T),
      dataset_details
    ),
    PORI_profile = c(
      "fixed_tree" = "PORI",
      "tree" = grep("constraint_tree_1",
                    grep(alignment_details$matrix_name,
                         grep(alignment_details$dataset, profile_topologies, value = T),
                         value = T),
                    value = T),
      dataset_details
    )
  )
  profile_commands <- lapply(1:length(PB_profiles),
                             create.phylobayes.profile.commands,
                             PB_profiles)
}


create.phylobayes.profile.commands <- function(j, PB_profiles){
  ## Assemble PhyloBayes commands
  profile_j <- PB_profiles[[j]]
  # pb_mpi -s matrix_40sp.fasta -T sister_to_lobopodia.tree -cat -poisson sister_to_lobopodia_chain1 (and chain_2)
  filename_prefix <- paste0(
    profile_j[["dataset"]],
    ".",
    profile_j[["matrix_name"]],
    ".PMSF_Phylobayes_CAT-POISSON.",
    profile_j[["fixed_tree"]]
  )
  pb_mpi_command <- paste0(
    profile_j[["pb_mpi_dir"]],
    "pb_mpi",
    " -s ",
    profile_j[["alignment"]],
    " -T ",
    profile_j[["tree"]],
    " -cat -poisson ",
    filename_prefix)
  pb_mpi_command_chain1 <- paste0(pb_mpi_command, ".chain1")
  pb_mpi_command_chain2 <- paste0(pb_mpi_command, ".chain2")
  # tracecomp -x 500 sister_to_lobopodia_chain1 sister_to_lobopodia_chain2
  tracecomp_command <- paste0(
    profile_j[["pb_mpi_dir"]],
    "tracecomp",
    " -x 500 ",
    filename_prefix, ".chain1",
    " ",
    filename_prefix, ".chain2"
  )
  # readpb_mpi -ss -x 500 10 sister_to_lobopodia_chain1
  readpb_mpi_command <- paste0(
    profile_j[["pb_mpi_dir"]],
    "readpb_mpi",
    " -ss -x 500 10 ",
    filename_prefix
  )
  readpb_mpi_command_chain1 <- paste0(readpb_mpi_command, ".chain1")
  readpb_mpi_command_chain2 <- paste0(readpb_mpi_command, ".chain2")
  # Output vector
  op_vector <- c(
    "dataset" = profile_j[["dataset"]],
    "matrix_name" = profile_j[["matrix_name"]],
    "alignment_path" = profile_j[["alignment"]],
    "topology" = profile_j[["fixed_tree"]],
    "ID" = filename_prefix,
    "pb_chain1" = pb_mpi_command_chain1,
    "pb_chain2" = pb_mpi_command_chain2,
    "tracecomp" = tracecomp_command,
    "readpb_chain1" = readpb_mpi_command_chain1,
    "readpb_chain2" = readpb_mpi_command_chain2
    )
  return(op_vector)
}

