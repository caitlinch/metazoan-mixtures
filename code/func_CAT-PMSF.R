## caitlinch/metazoan-mixtures/code/func_CAT-PMSF.R
# Functions for estimating CAT-PMSF models using PhyloBayes, and inferring trees with those models
# Caitlin Cherryh 2025



infer.phylobayes.profile.wrapper <- function(i, alignment_df, params, profile_topologies){
  # Return command lines to infer compositional profiles in PhyloBayes for a
  # given dataset

  i_alignment_df <- alignment_df[i, ]

  # Need to generate 2 profiles (under CTEN and PORI fixed topologies) with 2 chains each:
  i_dataset_details <- c(
    "dataset" = i_alignment_df$dataset,
    "matrix_name" = i_alignment_df$matrix_name,
    "alignment_dir" = params$alignments,
    "alignment_file" = i_alignment_df$alignment_path,
    "constraint_trees_dir" = params$constrained_CXX_trees,
    "phylobayes_dir" = params$phylobayes,
    "threads" = params$phylobayes_threads,
    "siteprofile_output_dir" = paste0(params$output_profiles, i_alignment_df$dataset, "/")
  )
  PB_profiles <- list(
    CTEN_profile = c(
      "fixed_tree" = "CTEN",
      "tree_id" = "ML_H1",
      "tree" = grep("ML_H1",
                    grep(i_dataset_details[["matrix_name"]],
                         grep(i_dataset_details[["dataset"]], profile_topologies, value = T),
                         value = T),
                    value = T),
      i_dataset_details
    ),
    PORI_profile = c(
      "fixed_tree" = "PORI",
      "tree_id" = "ML_H2",
      "tree" = grep("ML_H2",
                    grep(i_dataset_details[["matrix_name"]],
                         grep(i_dataset_details[["dataset"]], profile_topologies, value = T),
                         value = T),
                    value = T),
      i_dataset_details
    )
  )
  profile_command_df <- as.data.frame(
    do.call(
      rbind,
      lapply(
        1:length(PB_profiles),
        create.phylobayes.profile.commands,
        PB_profiles
      )
    )
  )
  return(profile_command_df)
}



create.phylobayes.profile.commands <- function(j, PB_profiles){
  ## Assemble PhyloBayes commands using file paths
  profile_j <- PB_profiles[[j]]
  # Create new output dir for each dataset
  mkdir_command <- paste0(
    "mkdir -p ",
    profile_j[["siteprofile_output_dir"]]
  )
  # From Giacomelli et al 2025:
  #     pb_mpi -s matrix_40sp.fasta -T sister_to_lobopodia.tree -cat -poisson sister_to_lobopodia_chain1
  # mpi command:
  #    mpirun -np {num threads} pb_mpi -s {fasta} -T {tree} -cat -poisson {chain_name}
  filename_prefix <- paste0(
    profile_j[["dataset"]],
    ".",
    profile_j[["matrix_name"]],
    ".PMSF_Phylobayes_CAT-POISSON.",
    profile_j[["fixed_tree"]]
  )
  # pb_mpi_command <- paste0(
  #   "mpirun",
  #   " -np ",
  #   profile_j[["threads"]],
  #   " pb_mpi -s ",
  #   profile_j[["alignment_dir"]],
  #   profile_j[["alignment_file"]],
  #   " -T ",
  #   profile_j[["constraint_trees_dir"]],
  #   profile_j[["tree"]],
  #   " -cat -poisson ",
  #   filename_prefix)
  pb_command <- paste0(
    profile_j[["phylobayes_dir"]], "pb",
    " -s ",
    profile_j[["alignment_dir"]],
    profile_j[["alignment_file"]],
    " -T ",
    profile_j[["constraint_trees_dir"]],
    profile_j[["tree"]],
    " -cat -poisson ",
    filename_prefix)
  pb_command_chain1 <- paste0(pb_command, ".chain1")
  pb_command_chain2 <- paste0(pb_command, ".chain2")
  # Giacomelli et al. 2025: tracecomp -x 500 sister_to_lobopodia_chain1 sister_to_lobopodia_chain2
  # mpi command:  mpirun -np {num threads} tracecomp -x 500 {chain1_name} {chain2_name}
  # tracecomp_mpi_command <- paste0(
  #   "mpirun",
  #   " -np ",
  #   profile_j[["threads"]],
  #   " tracecomp -x 500 ",
  #   profile_j[["siteprofile_output_dir"]],
  #   filename_prefix, ".chain1",
  #   " ",
  #   profile_j[["siteprofile_output_dir"]],
  #   filename_prefix, ".chain2"
  # )
  tracecomp_command <- paste0(
    profile_j[["phylobayes_dir"]], "tracecomp",
    " -x 500 ",
    profile_j[["siteprofile_output_dir"]],
    filename_prefix, ".chain1",
    " ",
    profile_j[["siteprofile_output_dir"]],
    filename_prefix, ".chain2"
  )
  mv_command <- paste0(
    "mv ",
    profile_j[["siteprofile_output_dir"]], "tracecomp.contdiff",
    " ",
    profile_j[["siteprofile_output_dir"]], filename_prefix, ".tracecomp.contdiff"
  )
  # Giacomelli et al. 2025: readpb_mpi -ss -x 500 10 sister_to_lobopodia_chain1
  # mpi command:  mpirun -np {num threads} readpb_mpi -ss -x 500 10 {chain1_name}
  # readpb_mpi_command <- paste0(
  #   "mpirun",
  #   " -np ", profile_j[["threads"]],
  #   " readpb_mpi -ss -x 500 10 ",
  #   profile_j[["siteprofile_output_dir"]],
  #   filename_prefix
  # )
  readpb_command <- paste0(
    profile_j[["phylobayes_dir"]], "readpb",
    " -ss -x 500 10 ",
    profile_j[["siteprofile_output_dir"]],
    filename_prefix
  )
  readpb_command_chain1 <- paste0(readpb_command, ".chain1")
  readpb_command_chain2 <- paste0(readpb_command, ".chain2")
  # Output vector
  op_vector <- c(
    "dataset" = profile_j[["dataset"]],
    "matrix_name" = profile_j[["matrix_name"]],
    "alignment_path" = profile_j[["alignment_file"]],
    "topology" = profile_j[["fixed_tree"]],
    "constraint_tree" = profile_j[["tree"]],
    "ID" = filename_prefix,
    "siteprofile_output_dir" = profile_j[["siteprofile_output_dir"]],
    "mkdir" = mkdir_command,
    "pb_chain1" = pb_command_chain1,
    "pb_chain2" = pb_command_chain2,
    "tracecomp" = tracecomp_command,
    "tracecomp_rename" = mv_command,
    "readpb_chain1" = readpb_command_chain1,
    "readpb_chain2" = readpb_command_chain2,
    "pb_siteprofile_file_chain1" = paste0(filename_prefix, ".chain1.siteprofiles"),
    "pb_siteprofile_file_chain2" = paste0(filename_prefix, ".chain2.siteprofiles")
  )
  return(op_vector)
}



convert.profile.for.iqtree <- function(i, commands_df, params){
  ## Generate command to convert PhyloBayes site distributions to IQ-TREE site frequencies
  commands_row <- commands_df[i, ]
  # convert-site-dists.py sister_to_lobopodia_chain1.siteprofiles
  convert_chain1 <- paste0(
    "python3 ",
    params$convert_site_dists,
    " ",
    commands_row$siteprofile_output_dir,
    commands_row$pb_siteprofile_file_chain1
  )
  convert_chain2 <- paste0(
    "python3 ",
    params$convert_site_dists,
    " ",
    commands_row$siteprofile_output_dir,
    commands_row$pb_siteprofile_file_chain2
  )
  sitefreqs_chain1 <- paste0(commands_row$ID, ".chain1.sitefreq")
  sitefreqs_chain2 <- paste0(commands_row$ID, ".chain2.sitefreq")
  op_vector <- c(as.character(commands_row),
                 convert_chain1, convert_chain2,
                 sitefreqs_chain1, sitefreqs_chain2)
  names(op_vector) <- c(names(commands_row),
                        "convert_siteprofiles_chain1", "convert_siteprofiles_chain2",
                        "sitefreqs_chain1", "sitefreqs_chain2")
  return(op_vector)
}


createPhyloBayesSlurmFile <- function(i, commands_df, jobscript_output_dir, jobscript_template){
  ## Convert row from dataframe into a runnable SLURM file for Dayhoff
  dataset_row <- commands_df[i, ]
  slurm_template <- readLines(jobscript_template)
  # Replace name of job with Dataset name
  dataset_slurm <- gsub(
    "--job-name=phylobayes_profiles",
    paste0("--job-name=pb_", dataset_row$dataset),
    slurm_template
  )
  # Add PhyloBayes commands
  dataset_slurm <- c(
    dataset_slurm,
    "# Make dataset directory",
    dataset_row$mkdir,
    "",
    "# Run PhyloBayes pb",
    dataset_row$pb_chain1,
    dataset_row$pb_chain2,
    "",
    "# Run PhyloBayes tracecomp",
    dataset_row$tracecomp,
    "",
    "# Run PhyloBayes readpb",
    dataset_row$readpb_chain1,
    dataset_row$readpb_chain2,
    ""
  )
  # Create output file path
  dataset_slurm_path <- paste0(
    jobscript_output_dir,
    dataset_row$ID,
    ".sh"
  )
  # Write file
  write(
    dataset_slurm,
    file = dataset_slurm_path
  )
  return(dataset_slurm_path)
}


createPythonConversionSlurmFile <- function(commands_df, jobscript_output_dir, jobscript_template){
  ## Create runnable SLURM file for use on Dayhoff to run the convert-site-dists.py
  #     script to transform .siteprofiles (PhyloBayes) files into .sitefrequencies
  #     (IQ-Tree) files
  slurm_template <- readLines(jobscript_template)
  ## Convert Chain 1 files
  slurm_chain1 <- c(
    gsub(
      "--job-name=phylobayes_profiles",
      "--job-name=sitefreqs_chain1",
      slurm_template
    ),
    "",
    paste0(
      "cd ",
      unique(dirname(commands_df$siteprofile_output_dir))
    ),
    "",
    commands_df$convert_siteprofiles_chain1
  )
  slurm_chain1_path <- paste0(
    jobscript_output_dir,
    "convert-site-dists.chain1.sh"
  )
  write(
    slurm_chain1,
    file = slurm_chain1_path
  )
  ## Convert chain 2 files
  slurm_chain2 <- c(
    gsub(
      "--job-name=phylobayes_profiles",
      "--job-name=sitefreqs_chain2",
      slurm_template
    ),
    "",
    paste0(
      "cd ",
      unique(dirname(commands_df$siteprofile_output_dir))
    ),
    "",
    commands_df$convert_siteprofiles_chain2
  )
  slurm_chain2_path <- paste0(
    jobscript_output_dir,
    "convert-site-dists.chain2.sh"
  )
  write(
    slurm_chain2,
    file = slurm_chain2_path
  )
  # Output
  op_files <- c(slurm_chain1_path, slurm_chain2_path)
  return(op_files)
}

