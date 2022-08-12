## caitlinch/metazoan-mixtures/func_constraint_trees.R
# Caitlin Cherryh 2022

# Functions for testing and applying constraint trees in iqtree2



estimate.ml.iqtree <- function(){
  # Function to call iqtree and estimate a maximum likelihood tree using best practices
}



run.iqtree.with.constraint.tree <- function(alignment_path, constraint_tree_file, partitioned_check = FALSE, partition_file = NA, 
                                            iqtree_path = "iqtree2", prefix = NA, model = NA, num_threads = 1){
  # Function to apply IQ-Tree to a series of alignments with a constraint tree
  
  # Set model for IQ-Tree run
  if (is.na(model) == TRUE){
    # If no model specified for IQ-Tree, use model finder (-m MFP) command
    iq_model = "MFP"
  } else {
    # Otherwise, use model specified in function call
    iq_model = model
  }
  
  # Add partition file if present
  if (partitioned_check == FALSE){
    partition_call <- ""
  } else if (partitioned_check == TRUE){
    partition_call <- paste0(" -p ", partition_file, " ")
  } 
  
  # Add prefix if present
  if (is.na(prefix) == TRUE){
    # If prefix is NA, do nothing
    prefix_call <- ""
  } else if (is.na(prefix) == FALSE){
    # If prefix is NA, add prefix to command line 
    prefix_call <- paste0(" --prefix ", prefix, " ")
  } 
  
  iqtree_call <- paste0(iqtree_path, " -s ", alignment_path,  partition_call, " -m ", iq_model, " -g ", constraint_tree_file, " -T ", num_threads,  prefix_call)
  
  # Run IQ-Tree
  system(iqtree_call)
}



run.one.constraint.tree <- function(index, df){
  # Quick function to take in a dataframe, take relevant variables, and call the run.iqtree.with.constraint.tree function
  
  # Identify row
  row <- df[index, ]
  
  # Feed row information into function call
  run.iqtree.with.constraint.tree(alignment_path = row$alignment_path, constraint_tree_file = row$constraint_tree_paths, 
                                  partitioned_check = row$partitioned_check, partition_file = row$partition_file, 
                                  iqtree_path = row$iqtree_path, prefix = row$constraint_prefixes, model = row$model,
                                  num_threads = row$num_threads)
}



run.one.constraint.dataframe <- function(csv_file){
  # Quick function to take in a dataframe, and estimate hypothesis trees by feeding it row by row into the run.one.constraint.tree function
  
  # Open the dataframe
  df <- read.csv(csv_file)
  # Estimate an ML tree in IQ-Tree for each constraint tree
  lapply(1:nrow(df), run.one.constraint.tree, df)
}



run.tree.mixture.model <- function(alignment_file, hypothesis_tree_file, partition_file, use.partition = FALSE, 
                                   prefix, model, number_parallel_threads, iqtree2_tree_mixtures_implementation){
  # Function runs the IQ-Tree2 mixture of trees model implementation given a sequence alignment, a set of hypothesis trees, and details about the model.
  
  # Add model if present
  if (is.na(model) == TRUE){
    # If no model specified for IQ-Tree, use model finder (-m MFP) command
    if (is.na(partition_file) == TRUE){
      # If no partition file is present, use only the MFP command
      model_call = "MFP"
    } else if (is.na(partition_file) == FALSE){
      # If a partition file is present, use MFP+MERGE
      model_call = "MFP+MERGE"
    }
    # Extend the model to have the +TR command 
    model_call = paste0("'",model_call, "+TR'")
  } else if (is.na(model) == FALSE){
    # If a model is provided, use that model
    model_call = model
    # Extend the model to have the +TR command 
    model_call = paste0("'",model_call, "+TR'")
  }
  
  # Add partition file if present
  if (is.na(partition_file) == TRUE){
    # If partition_file is NA, do nothing
    partition_call <- ""
  } else if (is.na(partition_file) == FALSE){
    # If prefix is NA, add prefix to command line 
    partition_call <- paste0(" -p ", partition_file, " ")
  } 
  
  # Add prefix if present
  if (is.na(prefix) == TRUE){
    # If prefix is NA, do nothing
    prefix_call <- ""
  } else if (is.na(prefix) == FALSE){
    # If prefix is NA, add prefix to command line 
    prefix_call <- paste0(" -pre ", prefix, " ")
  }
  
  if (use.partition == FALSE){
    # Assemble the command for the tree mixtures model
    treemix_command <- paste0(iqtree2_tree_mixtures_implementation, " -s ", alignment_file, " -m  ", model_call, 
                              " -te ", hypothesis_tree_file, " -nt ", number_parallel_threads, prefix_call)
  } else if (use.partition == TRUE){
    # Assemble the command for the tree mixtures model
    treemix_command <- paste0(iqtree2_tree_mixtures_implementation, " -s ", alignment_file, partition_call, " -m ", model_call, 
                              " -te ", hypothesis_tree_file, " -nt ", number_parallel_threads, prefix_call)
  }
  
  # Change working directories (to store IQ-Tree output files in the right place)
  setwd(dirname(hypothesis_tree_file))
  
  # Call IQ-Tree2 with the command
  system(treemix_command)
  
}



create.constraint.trees <- function(dataset, dataset_constraint_tree_dir, model, model_id, outgroup_taxa, ctenophora_taxa, 
                                    porifera_taxa, sponges_1_taxa, sponges_2_taxa, placozoa_taxa, cnidaria_taxa, bilateria_taxa,
                                    alignment_file, partitioned_check, partition_file, iqtree_path){
  # Function to create the constraint trees and constraint tree information data frame, for a given dataset and model
  
  ## Hypothesis 1: Ctenophora-sister
  # Tree: (outgroup_taxa, (ctenophora_taxa, (porifera_taxa, (placozoa_taxa, cnidaria_taxa, bilateria_taxa))))
  # Construct constraint tree
  constraint_tree_1 <- paste0("((", 
                              paste(outgroup_taxa, collapse = ", "), 
                              "),((", 
                              paste(ctenophora_taxa, collapse = ", "), 
                              "),((", 
                              paste(c(porifera_taxa), collapse = ", "),
                              "), (",
                              paste(c(placozoa_taxa, cnidaria_taxa, bilateria_taxa), collapse = ", "),
                              "))));")
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, dataset, "_", model_id,"_constraint_tree_", "1", ".nex")
  write(constraint_tree_1, file = constraint_tree_file_name)
  
  
  ## Hypothesis 2: Porifera-sister
  # Tree: (outgroup_taxa, (porifera_taxa, (ctenophora_taxa, (placozoa_taxa, cnidaria_taxa, bilateria_taxa))))
  # Construct constraint tree
  constraint_tree_2 <- paste0("((", 
                              paste(outgroup_taxa, collapse = ", "), 
                              "),((", 
                              paste(porifera_taxa, collapse = ", "), 
                              "),((", 
                              paste(c(ctenophora_taxa), collapse = ", "),
                              "),(",
                              paste(c(placozoa_taxa, cnidaria_taxa, bilateria_taxa), collapse = ", "), 
                              "))));")
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, dataset, "_", model_id,"_constraint_tree_", "2", ".nex")
  write(constraint_tree_2, file = constraint_tree_file_name)
  
  
  ## Hypothesis 3: Porifera+Ctenophora-sister
  # Tree: (outgroup_taxa, ((porifera_taxa, ctenophora_taxa), (placozoa_taxa, (cnidaria_taxa, bilateria_taxa))))
  # Construct constraint tree
  constraint_tree_3 <- paste0("((", 
                              paste(outgroup_taxa, collapse = ", "), 
                              "),((", 
                              paste(c(porifera_taxa, ctenophora_taxa), collapse = ", "), 
                              "),(", 
                              paste(c(placozoa_taxa), collapse = ", "),
                              ", (",
                              paste(c(cnidaria_taxa, bilateria_taxa), collapse = ", "), 
                              "))));")
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, dataset, "_", model_id,"_constraint_tree_", "3", ".nex")
  write(constraint_tree_3, file = constraint_tree_file_name)
  
  ## Hypothesis 4: Paraphyletic sponges, Porifera-sister
  # Tree: (outgroup_taxa, (sponges_1_taxa, (sponges_2_taxa, (ctenophora_taxa, (placozoa_taxa, cnidaria_taxa, bilateria_taxa)))))
  # Construct constraint tree
  constraint_tree_4 <- paste0("((", 
                              paste(outgroup_taxa, collapse = ", ") ,
                              ") ,((", 
                              paste(sponges_1_taxa, collapse = ", "), 
                              "), ((", 
                              paste(sponges_2_taxa, collapse = ", "), 
                              "), ((", 
                              paste(c(ctenophora_taxa), collapse = ", "),
                              "), (",
                              paste(c(placozoa_taxa, cnidaria_taxa, bilateria_taxa), collapse = ", "),
                              ")))));")
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, dataset, "_", model_id,"_constraint_tree_", "4", ".nex")
  write(constraint_tree_4, file = constraint_tree_file_name)
  
  ## Hypothesis 5: Paraphyletic sponges, Ctenophora-sister
  # Tree: (outgroup_taxa, (ctenophora_taxa, (sponges_1_taxa, (sponges_2_taxa, (placozoa_taxa, cnidaria_taxa, bilateria_taxa)))))
  # Construct constraint tree
  constraint_tree_5 <- paste0("((", 
                              paste(outgroup_taxa, collapse = ", "),
                              ") ,((",
                              paste(ctenophora_taxa, collapse = ", "),
                              "), ((", 
                              paste(sponges_1_taxa, collapse = ", "),
                              "), ((", 
                              paste(c(sponges_2_taxa), collapse = ", "),
                              "), (",
                              paste(c(placozoa_taxa, cnidaria_taxa, bilateria_taxa), collapse = ", "),
                              ")))));")
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, dataset, "_", model_id,"_constraint_tree_", "5", ".nex")
  write(constraint_tree_5, file = constraint_tree_file_name)
  
  ## Hypothesis 6: Paraphyletic sponges, Porifera-sister
  # Tree: (outgroup_taxa, (sponges_2_taxa, (sponges_1_taxa, (ctenophora_taxa, (placozoa_taxa, cnidaria_taxa, bilateria_taxa)))))
  # Construct constraint tree
  constraint_tree_6 <- paste0("((", 
                              paste(outgroup_taxa, collapse = ", ") ,
                              ") ,((", 
                              paste(sponges_2_taxa, collapse = ", "), 
                              "), ((", 
                              paste(sponges_1_taxa, collapse = ", "), 
                              "), ((", 
                              paste(c(ctenophora_taxa), collapse = ", "),
                              "), (",
                              paste(c(placozoa_taxa, cnidaria_taxa, bilateria_taxa), collapse = ", "),
                              ")))));")
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, dataset, "_", model_id,"_constraint_tree_", "6", ".nex")
  write(constraint_tree_6, file = constraint_tree_file_name)
  
  ## Hypothesis 7: Paraphyletic sponges, Ctenophora-sister
  # Tree: (outgroup_taxa, (ctenophora_taxa, (sponges_2_taxa, (sponges_1_taxa, (placozoa_taxa, cnidaria_taxa, bilateria_taxa)))))
  # Construct constraint tree
  constraint_tree_7 <- paste0("((", 
                              paste(outgroup_taxa, collapse = ", "),
                              ") ,((",
                              paste(ctenophora_taxa, collapse = ", "),
                              "), ((", 
                              paste(sponges_2_taxa, collapse = ", "),
                              "), ((", 
                              paste(c(sponges_1_taxa), collapse = ", "),
                              "), (",
                              paste(c(placozoa_taxa, cnidaria_taxa, bilateria_taxa), collapse = ", "),
                              ")))));")
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, dataset, "_", model_id,"_constraint_tree_", "7", ".nex")
  write(constraint_tree_7, file = constraint_tree_file_name)
  
  # Assemble dataframe of information about the constraint trees
  constraint_df <- data.frame(constraint_tree_id = 1:7,
                              constraint_tree_paths = paste0(dataset_constraint_tree_dir, dataset, "_", model_id,"_constraint_tree_", 1:7, ".nex"),
                              constraint_prefixes = paste0(dataset, "_", model, "_ML_H", 1:7),
                              alignment_path = alignment_file,
                              model = model,
                              iqtree_path = iqtree_path,
                              constraint_trees = c(constraint_tree_1, constraint_tree_2, constraint_tree_3, 
                                                   constraint_tree_4, constraint_tree_5, constraint_tree_6,
                                                   constraint_tree_7),
                              num_threads = number_parallel_threads,
                              partitioned = partitioned_check,
                              partition_file = partition_file)
  
  # Write dataframe of information about constraint trees
  constraint_df_path <- paste0(dataset_constraint_tree_dir, dataset, "_", model_id, "_constraint_tree_parameters.csv")
  write.csv(constraint_df, constraint_df_path, row.names = FALSE)
  
  # Return the constraint tree dataframe
  return(constraint_df)
}
