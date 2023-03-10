## caitlinch/metazoan-mixtures/func_constraint_trees.R
# Caitlin Cherryh 2022

# Functions for estimating maximum likelihood trees and constrained maximum likelihood trees,
#     and for running the Mixture of Trees and Sites (MAST) model in iqtree2





#### Estimating an ML tree in IQ-Tree with a specified model ####
estimate.ml.iqtree <- function(iqtree2, alignment_file, model = "MFP", mset = NA, mrate = NA, partition_file = NA, 
                               prefix = NA, number_parallel_threads = "AUTO", number_of_bootstraps = NA,
                               redo = FALSE, safe = FALSE, run.iqtree = TRUE){
  # Function to call iqtree and estimate a maximum likelihood tree using best practices
  
  # Add partition file if present
  if (is.na(partition_file) == TRUE){
    # If the partition file is NA, there is no partition file for this alignment
    partition_call <- ""
    # Check whether model is specified
    if (is.na(model) == TRUE){
      # Tell IQ-Tree to use ModelFinder
      model_call = ""
    } else if (is.na(model) == FALSE){
      # Do not use ModelFinder
      model_call = paste0(" -m ", model, " ")
    }
    # Check whether mset is specified
    if (is.na(mset) == TRUE){
      # If mset = NA, then no mset option is specified.
      mset_call = ""
      # Tell IQ-Tree to use ModelFinder
    } else if (is.na(mset) == FALSE){
      # If mset is specified, add mset command
      mset_call <- paste0(" -mset '", mset, "' ")
    }
    # Check whether mrate is specified
    if (is.na(mrate) == TRUE){
      # If mrate = NA, then no mset option is specified.
      mrate_call = ""
      # Tell IQ-Tree to use ModelFinder
    } else if (is.na(mrate) == FALSE){
      # If mrate is specified, add mrate command
      mrate_call <- paste0(" -mrate '", mrate, "' ")
    }
  } else if (is.na(partition_file) == FALSE){
    # If the partition file is not NA, add the command for a partition file to the command line for iqtree
    partition_call <- paste0(" -p ", partition_file, " ")
    # If there is a partition file, set the model selection to include a merging step
    model_call = " -m MFP+MERGE "
    # There is no mset command (models are already specified in the partition file)
    mset_call = ""
    mrate_call = ""
  } 
  
  # If prefix is specified, add a prefix command to the command line
  if (is.na(prefix) == FALSE){
    prefix_call = paste0(" -pre ", prefix, " ")
  } else if (is.na(prefix) == TRUE){
    prefix_call = ""
  }
  
  # If number of bootstraps is specified, add a bootstrap command to the command line
  if (is.na(number_of_bootstraps) == FALSE){
    bootstrap_call = paste0(" -bb ", number_of_bootstraps, " ")
  } else if (is.na(number_of_bootstraps) == TRUE){
    bootstrap_call = ""
  }
  
  # If redo is TRUE, add redo command to IQ-Tree call
  if (redo == TRUE){
    redo_command = " -redo "
  } else if (redo == FALSE){
    redo_command = ""
  }
  
  # If safe is TRUE, add safe command to IQ-Tree call
  if (safe == TRUE){
    safe_command = " -safe "
  } else if (safe == FALSE){
    safe_command = ""
  }
  
  # Assemble the command line
  iqtree_call <- paste0(iqtree2, " -s ", alignment_file, partition_call, model_call, mset_call, mrate_call,
                        bootstrap_call, redo_command, safe_command, " -nt ", number_parallel_threads,
                        prefix_call)
  
  if (run.iqtree == TRUE){
    # Call iqtree to estimate the tree
    system(iqtree_call)
  } # end if (run.iqtree == TRUE)
  
  # Return the iqtree2 command
  return(iqtree_call)
} # end function


ml.iqtree.wrapper <- function(i, iqtree_path, df){
  # Function to take row from a dataframe and call estimate.ml.iqtree using information from that row
  
  # Extract the row
  row <- df[i,]
  
  # Create the command line for iqtree
  iqtree_call <- estimate.ml.iqtree(iqtree2 = iqtree_path, alignment_file = row$alignment_file, 
                                    model = row$model_m, mset = row$model_mset, mrate = row$model_mrate, partition_file = NA, 
                                    prefix = row$prefix, number_parallel_threads = row$iqtree_num_threads, 
                                    number_of_bootstraps = row$iqtree_num_bootstraps,
                                    redo = FALSE, safe = FALSE, run.iqtree = FALSE)
  
  # Return the command line for iqtree
  return(iqtree_call)
} # end function





#### Prepare to estimate constraint trees ####
determine.best.ML.model.wrapper <- function(row_id, completed_runs_df, ML_output_df){
  # Function to take in a row number, the dataframe of which runs have completed, and the data frame of output from the
  #     maximum likelihood runs, and feed that information into another function to pick the best ML models
  
  # Find the row
  row <- completed_runs_df[row_id, ]
  # Reduce the output dataframe to just the dataset/matrix in that row
  dataset_output_df <- ML_output_df[ML_output_df$dataset == row$dataset & ML_output_df$matrix_name == row$matrix_name,]
  # Input information from the row into the determine.best.ML.model function
  dataset_best_model_df <- determine.best.ML.model(dataset = row$dataset, matrix = row$matrix_name, dataset_output_df = dataset_output_df)
  # Return the best model(s) for the dataset
  return(dataset_best_model_df)
} # end function


determine.best.ML.model <- function(dataset, matrix, dataset_output_df){
  # Function to determine the best substitution model (by BIC) for a particular dataset and alignment combination,
  #     using the output from the maximum likelihood runs
  
  # Sort the dataset by BIC (lowest to highest - remember lower is better)
  dataset_output_df <- dataset_output_df[order(dataset_output_df$best_model_BIC, decreasing = FALSE),]
  # Determine which model has the lowest BIC
  model_lowest_BIC <- dataset_output_df$model_code[1]
  # Check whether ModelFinder has the lowest BIC
  is.ModelFinder.BIC.best <- (model_lowest_BIC == "ModelFinder")
  # Prepare the output
  if (is.ModelFinder.BIC.best == TRUE){
    # If ModelFinder has the best BIC, return the ModelFinder row only
    best_model_df <- dataset_output_df[which(dataset_output_df$model_code == "ModelFinder"),]
  } else if (is.ModelFinder.BIC.best == FALSE){
    # If ModelFinder does not have the best BIC, return both the ModelFinder row and the row for the model with the best BIC
    best_model_df <- dataset_output_df[c(1, which(dataset_output_df$model_code == "ModelFinder")),]
  }
  # Return the dataframe of information from the best model(s)
  return(best_model_df)
} # end function





#### Creating constraint trees ####
constraint.tree.wrapper <- function(i, output_directory, iqtree_path, iqtree_num_threads = 1, 
                                    dataset_info = all_datasets, matrix_taxa_info = matrix_taxa, 
                                    ml_output_df, ml_tree_tips_df, force.update.constraint.trees = TRUE){
  
  # Function to take a row from the ML output dataframe and create the constraint trees plus parameters to estimate the hypothesis trees
  
  ## Extract the relevant row
  row <- ml_output_df[i, ]
  
  ## Extract the relevant list of taxa for this dataframe
  # First, check whether this matrix is included in the keys of the matrix_taxa list
  row_dataset <- row$dataset
  row_key <- paste0(row$dataset, ".", row$matrix_name, ".", row$sequence_format)
  list_keys <- names(matrix_taxa_info)
  # Check if row_key in list_key
  if ((row_key %in% list_keys) == FALSE){
    # If row key is not in list key, then all taxa for this dataset have the same names
    # Extract the object containing those taxa names
    constraint_clades <- dataset_info[[row_dataset]]
  } else if ((row_key %in% list_keys) == TRUE){
    # First, identify the list of taxa in this matrix
    keep_taxa <- matrix_taxa_info[[row_key]]
    # Secondly, extract the taxa clades for this dataset
    dataset_taxa_clades <- dataset_info[[row_dataset]]
    # Make a copy of the clades object
    constraint_clades <- dataset_taxa_clades
    # Lastly, remove any taxa that is NOT in the keep_taxa from the constraint clades
    #   i.e., remove any taxa from this dataset that are NOT present in this matrix
    #   (as some datasets have multiple matrices, with different taxon sampling or different taxon naming conventions)
    constraint_clades$Bilateria <- dataset_taxa_clades$Bilateria[which(dataset_taxa_clades$Bilateria %in% keep_taxa)]
    constraint_clades$Cnidaria <- dataset_taxa_clades$Cnidaria[which(dataset_taxa_clades$Cnidaria %in% keep_taxa)]
    constraint_clades$Placozoa <- dataset_taxa_clades$Placozoa[which(dataset_taxa_clades$Placozoa %in% keep_taxa)]
    constraint_clades$Porifera <- dataset_taxa_clades$Porifera[which(dataset_taxa_clades$Porifera %in% keep_taxa)]
    constraint_clades$Ctenophora <- dataset_taxa_clades$Ctenophora[which(dataset_taxa_clades$Ctenophora %in% keep_taxa)]
    constraint_clades$Outgroup <- dataset_taxa_clades$Outgroup[which(dataset_taxa_clades$Outgroup %in% keep_taxa)]
    constraint_clades$Sponges_Calcarea <- dataset_taxa_clades$Sponges_Calcarea[which(dataset_taxa_clades$Sponges_Calcarea %in% keep_taxa)]
    constraint_clades$Sponges_Homoscleromorpha <- dataset_taxa_clades$Sponges_Homoscleromorpha[which(dataset_taxa_clades$Sponges_Homoscleromorpha %in% keep_taxa)]
    constraint_clades$Sponges_Hexactinellida <- dataset_taxa_clades$Sponges_Hexactinellida[which(dataset_taxa_clades$Sponges_Hexactinellida %in% keep_taxa)]
    constraint_clades$Sponges_Demospongiae <- dataset_taxa_clades$Sponges_Demospongiae[which(dataset_taxa_clades$Sponges_Demospongiae %in% keep_taxa)]
  }
  
  ## Remove any taxa from the constraint_clades that are not included in the ML tree for the alignment
  # Create the unique identifier for this row
  unique_id <- paste0(row$dataset, ".", row$matrix_name)
  # Extract the column of tip labels from the relevant unique_id column of the ml_tree_tips_df
  tree_tips_raw <- ml_tree_tips_df[[c(unique_id)]]
  tree_tips_cleaned <- na.omit(tree_tips_raw)
  tree_tips <- as.character(tree_tips_cleaned)
  # Check each of the clades and remove any tips not in the list of tree tips
  constraint_clades$Bilateria <- constraint_clades$Bilateria[(constraint_clades$Bilateria %in% tree_tips)]
  constraint_clades$Cnidaria <- constraint_clades$Cnidaria[(constraint_clades$Cnidaria %in% tree_tips)]
  constraint_clades$Placozoa <- constraint_clades$Placozoa[(constraint_clades$Placozoa %in% tree_tips)]
  constraint_clades$Porifera <- constraint_clades$Porifera[(constraint_clades$Porifera %in% tree_tips)]
  constraint_clades$Ctenophora <- constraint_clades$Ctenophora[(constraint_clades$Ctenophora %in% tree_tips)]
  constraint_clades$Outgroup <- constraint_clades$Outgroup[(constraint_clades$Outgroup %in% tree_tips)]
  constraint_clades$Sponges_Calcarea <- constraint_clades$Sponges_Calcarea[(constraint_clades$Sponges_Calcarea %in% tree_tips)]
  constraint_clades$Sponges_Homoscleromorpha <- constraint_clades$Sponges_Homoscleromorpha[(constraint_clades$Sponges_Homoscleromorpha %in% tree_tips)]
  constraint_clades$Sponges_Hexactinellida <- constraint_clades$Sponges_Hexactinellida[(constraint_clades$Sponges_Hexactinellida %in% tree_tips)]
  constraint_clades$Sponges_Demospongiae <- constraint_clades$Sponges_Demospongiae[(constraint_clades$Sponges_Demospongiae %in% tree_tips)]
  
  ## Create the constraint tree dataframe by calling the create.constraint.trees function
  constraint_df <- create.constraint.trees(dataset = row$dataset, 
                                           matrix_name = row$matrix_name,
                                           model_code = row$model_code,
                                           prefix = row$prefix, 
                                           dataset_constraint_tree_dir = output_directory,
                                           best_model = row$best_model,
                                           estimated_rates = row$estimated_rates,
                                           estimated_gamma = row$estimated_gamma,
                                           estimated_state_frequencies = row$estimated_state_frequencies,
                                           constraint_clades = constraint_clades,
                                           alignment_file = row$alignment_file, 
                                           partitioned_check = FALSE, 
                                           partition_file = NA, 
                                           iqtree_path = iqtree_path, 
                                           number_parallel_threads = iqtree_num_threads,
                                           force.update.constraint.trees = TRUE)
  ## Return the constraint tree dataframe
  return(constraint_df)
} # end function


format.constraint.tree.clade <- function(clade){
  # Function to take in a vector of species and return a nicely formatted character object to paste into a constraint tree
  
  # Check how many taxa are in the clade
  clade_size <- length(clade)
  
  # Format the clade
  if (clade_size == 1){
    clade_formatted = clade
  } else if (clade_size > 1){
    clade_formatted = paste0("(", paste(clade, collapse = ", "), ")")
  } else if (clade_size < 1){
    clade_formatted <- ""
  }
  
  # Return the nicely formatted clade
  return(clade_formatted)
}


create.constraint.trees <- function(dataset, matrix_name, model_code, prefix = NA, dataset_constraint_tree_dir, 
                                    best_model, estimated_rates = NA, estimated_gamma = NA, estimated_state_frequencies = NA, 
                                    constraint_clades, alignment_file, partitioned_check = FALSE, partition_file = NA, 
                                    iqtree_path, number_parallel_threads = 1, force.update.constraint.trees = TRUE){
  # Function to create the constraint trees and constraint tree information data frame, for a given dataset and model
  # Does not include Placozoa - Placozoa is not relevant to the question I am investigating (and is only ~1 taxa in most data sets)
  
  ### Prepare to construct constraint trees
  # Make sure you have an output id, which is a unique identifier for each dataset/alignment/model combination.
  if (is.na(prefix) == FALSE){
    # If a prefix is provided, use it in the file names
    prefix = prefix
  } else if (is.na(prefix) == TRUE){
    # If no prefix is provided, create one
    prefix <- paste0(dataset, ".", matrix_name, ".", model_code)
  }
  
  ### Construct constraint trees
  # Generate file names for all 5 constraint trees
  constraint_tree_1_file_name <- paste0(dataset_constraint_tree_dir, dataset, ".", matrix_name, ".constraint_tree_", "1", ".nex")
  constraint_tree_2_file_name <- paste0(dataset_constraint_tree_dir, dataset, ".", matrix_name, ".constraint_tree_", "2", ".nex")
  constraint_tree_3_file_name <- paste0(dataset_constraint_tree_dir, dataset, ".", matrix_name, ".constraint_tree_", "3", ".nex")
  constraint_tree_4_file_name <- paste0(dataset_constraint_tree_dir, dataset, ".", matrix_name, ".constraint_tree_", "4", ".nex")
  constraint_tree_5_file_name <- paste0(dataset_constraint_tree_dir, dataset, ".", matrix_name, ".constraint_tree_", "5", ".nex")
  
  # Split the taxa into clades
  outgroup_taxa = constraint_clades$Outgroup
  ctenophora_taxa = constraint_clades$Ctenophora
  porifera_taxa = constraint_clades$Porifera
  sponges_1_taxa = as.character(unlist(constraint_clades[c(constraint_clades$Sponges_1)]))
  sponges_2_taxa = as.character(unlist(constraint_clades[c(constraint_clades$Sponges_2)])) 
  cnidaria_taxa = constraint_clades$Cnidaria
  bilateria_taxa = constraint_clades$Bilateria
  
  # Format the outgroup clades nicely
  outgroup_taxa_formatted <- format.constraint.tree.clade(outgroup_taxa)
  ctenophora_taxa_formatted <- format.constraint.tree.clade(ctenophora_taxa)
  porifera_taxa_formatted <- format.constraint.tree.clade(porifera_taxa)
  ctenophora_porifera_taxa_formatted <- format.constraint.tree.clade(c(ctenophora_taxa, porifera_taxa))
  sponges_1_taxa_formatted <- format.constraint.tree.clade(sponges_1_taxa)
  sponges_2_taxa_formatted <- format.constraint.tree.clade(sponges_2_taxa)
  cnidaria_taxa_formatted <- format.constraint.tree.clade(cnidaria_taxa)
  bilateria_taxa_formatted <- format.constraint.tree.clade(bilateria_taxa)
  cnidaria_bilateria_taxa_formatted <- format.constraint.tree.clade(c(cnidaria_taxa, bilateria_taxa))
  
  
  # Only need to create one set of hypothesis trees per dataset/matrix combination - create constraint tree if it doesn't exist
  ### Hypothesis 1: Ctenophora-sister
  # Tree: (outgroup_taxa, (ctenophora_taxa, (porifera_taxa, (cnidaria_taxa, bilateria_taxa))));
  if ((file.exists(constraint_tree_1_file_name) == FALSE) | (force.update.constraint.trees == TRUE)){
    ## Construct constraint tree
    constraint_tree_1 <- paste0("(", outgroup_taxa_formatted, ", (", ctenophora_taxa_formatted, ", (", porifera_taxa_formatted, ", ", cnidaria_bilateria_taxa_formatted, ")));")
    write(constraint_tree_1, file = constraint_tree_1_file_name)
  } else if (file.exists(constraint_tree_1_file_name) == TRUE){ 
    constraint_tree_1 <- readLines(constraint_tree_1_file_name)
  }
  
  ### Hypothesis 2: Porifera-sister
  # Tree: (outgroup_taxa, (porifera_taxa, (ctenophora_taxa, (cnidaria_taxa, bilateria_taxa))));
  if ((file.exists(constraint_tree_2_file_name) == FALSE) | (force.update.constraint.trees == TRUE)){
    ## Construct constraint tree
    constraint_tree_2 <- paste0("(", outgroup_taxa_formatted, ", (", porifera_taxa_formatted, ", (", ctenophora_taxa_formatted, ", ", cnidaria_bilateria_taxa_formatted, ")));")
    write(constraint_tree_2, file = constraint_tree_2_file_name)
  } else if (file.exists(constraint_tree_2_file_name) == TRUE){ 
    constraint_tree_2 <- readLines(constraint_tree_2_file_name)
  }
  
  ### Hypothesis 3: (Ctenophore+Porifera)-sister
  # Tree: (outgroup_taxa, ((porifera_taxa, ctenophora_taxa), (cnidaria_taxa, bilateria_taxa)));
  if ((file.exists(constraint_tree_3_file_name) == FALSE) | (force.update.constraint.trees == TRUE)){
    ## Construct constraint tree
    constraint_tree_3 <- paste0("(", outgroup_taxa_formatted, ", (", ctenophora_porifera_taxa_formatted, ", ", cnidaria_bilateria_taxa_formatted, "));")
    write(constraint_tree_3, file = constraint_tree_3_file_name)
  } else if (file.exists(constraint_tree_3_file_name) == TRUE){ 
    constraint_tree_3 <- readLines(constraint_tree_3_file_name)
  }
  
  # Only create a constraint tree for Hypothesis 4 and 5 IF there is at least one taxa in each sponge group (sponges_1 and sponges_2)
  if (length(sponges_1_taxa) > 0 & length(sponges_2_taxa) > 0){
    # There is >= 1 taxa in each sponge group (sponges_1 and sponges_2). Estimate constraint trees 4 and 5.
    
    ### Hypothesis 4: Paraphyletic sponges, Ctenophora-sister
    # Tree: (outgroup_taxa, (ctenophora_taxa, (sponges_1_taxa, (sponges_2_taxa, (cnidaria_taxa, bilateria_taxa)))));
    if ((file.exists(constraint_tree_4_file_name) == FALSE) | (force.update.constraint.trees == TRUE)){
      ## Construct constraint tree
      constraint_tree_4 <- paste0("(", outgroup_taxa_formatted, ", (", ctenophora_taxa_formatted, ", (", sponges_1_taxa_formatted, ", (", sponges_2_taxa_formatted,  ", ", cnidaria_bilateria_taxa_formatted, "))));")
      write(constraint_tree_4, file = constraint_tree_4_file_name)
    } else if (file.exists(constraint_tree_4_file_name) == TRUE){ 
      constraint_tree_4 <- readLines(constraint_tree_4_file_name)
    }
    
    ### Hypothesis 5: Paraphyletic sponges, Porifera-sister
    # Tree: (outgroup_taxa, (sponges_1_taxa, (sponges_2_taxa, (ctenophora_taxa, (cnidaria_taxa, bilateria_taxa)))));
    if ((file.exists(constraint_tree_5_file_name) == FALSE) | (force.update.constraint.trees == TRUE)){
      ## Construct constraint tree
      constraint_tree_5 <- paste0("(", outgroup_taxa_formatted, ", (", sponges_1_taxa_formatted, ", (", sponges_2_taxa_formatted, ", (", ctenophora_taxa_formatted,  ", ", cnidaria_bilateria_taxa_formatted, "))));")
      write(constraint_tree_5, file = constraint_tree_5_file_name)
    } else if (file.exists(constraint_tree_5_file_name) == TRUE){ 
      constraint_tree_5 <- readLines(constraint_tree_5_file_name)
    }
    
    # Assemble dataframe of information about the 5 constraint trees
    constraint_df <- data.frame(dataset = dataset,
                                matrix_name = matrix_name,
                                model_code = model_code,
                                prefix = prefix,
                                best_model = best_model,
                                estimated_rates = estimated_rates,
                                estimated_gamma = estimated_gamma,
                                estimated_state_frequencies = estimated_state_frequencies,
                                constraint_tree_hypothesis = c("Ctenophora-sister", "Porifera-sister", "(Ctenophora+Porifera)-sister", "Ctenophora-sister (Paraphyletic sponges)", "Porifera-sister (Paraphyletic sponges)"), 
                                constraint_tree_number = 1:5,
                                constraint_tree_id = paste0(prefix, ".ML_C", 1:5),
                                constraint_tree_paths = c(constraint_tree_1_file_name, constraint_tree_2_file_name, constraint_tree_3_file_name, constraint_tree_4_file_name, constraint_tree_5_file_name),
                                hypothesis_tree_prefixes = paste0(prefix, ".ML_H", 1:5),
                                alignment_path = alignment_file,
                                iqtree_path = iqtree_path,
                                constraint_trees = c(constraint_tree_1, constraint_tree_2, constraint_tree_3, constraint_tree_4, constraint_tree_5),
                                num_threads = number_parallel_threads,
                                partitioned = partitioned_check,
                                partition_file = partition_file)
  } else {
    # There is not >= 1 taxa in each sponge group (sponges_1 and sponges_2), so skip constraint trees 4 and 5
    # Assemble dataframe of information about the 3 constraint trees
    constraint_df <- data.frame(dataset = dataset,
                                matrix_name = matrix_name,
                                model_code = model_code,
                                prefix = prefix,
                                best_model = best_model,
                                estimated_rates = estimated_rates,
                                estimated_gamma = estimated_gamma,
                                estimated_state_frequencies = estimated_state_frequencies,
                                constraint_tree_hypothesis = c("Ctenophora-sister", "Porifera-sister", "(Ctenophora+Porifera)-sister"), 
                                constraint_tree_number = 1:5,
                                constraint_tree_id = paste0(prefix, ".ML_C", 1:3),
                                constraint_tree_paths = c(constraint_tree_1_file_name, constraint_tree_2_file_name, constraint_tree_3_file_name),
                                hypothesis_tree_prefixes = paste0(prefix, ".ML_H", 1:3),
                                alignment_path = alignment_file,
                                iqtree_path = iqtree_path,
                                constraint_trees = c(constraint_tree_1, constraint_tree_2, constraint_tree_3),
                                num_threads = number_parallel_threads,
                                partitioned = partitioned_check,
                                partition_file = partition_file)
  }
  
  # Write dataframe of information about constraint trees
  constraint_df_path <- paste0(dataset_constraint_tree_dir, prefix, "_constraint_tree_parameters.csv")
  write.csv(constraint_df, constraint_df_path, row.names = FALSE)
  
  # Return the constraint tree dataframe
  return(constraint_df)
} # end function


create.constraint.trees.Placozoa <- function(dataset, prefix = NA, dataset_constraint_tree_dir, best_model, model_code, outgroup_taxa, ctenophora_taxa, 
                                             porifera_taxa, sponges_1_taxa, sponges_2_taxa, placozoa_taxa, cnidaria_taxa, bilateria_taxa,
                                             alignment_file, partitioned_check, partition_file, iqtree_path, number_parallel_threads){
  # Function to create the constraint trees and constraint tree information data frame, for a given dataset and model
  # Includes Placozoa
  
  # Make sure you have an output id, which is a unique identifier for each dataset/alignment/model combination.
  if (is.na(prefix) == FALSE){
    # If a prefix is provided, use it in the file names
    prefix = prefix
  }
  else if (is.na(prefix) == TRUE){
    # If no prefix is provided, create one
    prefix <- paste0(dataset, ".", model_code)
  }
  
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
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, prefix, "_constraint_tree_", "1", ".nex")
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
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, prefix, "_constraint_tree_", "2", ".nex")
  write(constraint_tree_2, file = constraint_tree_file_name)
  
  # Assemble dataframe of information about the constraint trees
  constraint_df <- data.frame(dataset = dataset,
                              prefix = prefix,
                              model_code = model_code,
                              constraint_tree_id = 1:2,
                              constraint_tree_paths = paste0(dataset_constraint_tree_dir, prefix, "_constraint_tree_", 1:2, ".nex"),
                              constraint_prefixes = paste0(prefix, "_ML_H", 1:2),
                              alignment_path = alignment_file,
                              best_model = best_model,
                              iqtree_path = iqtree_path,
                              constraint_trees = c(constraint_tree_1, constraint_tree_2),
                              num_threads = number_parallel_threads,
                              partitioned = partitioned_check,
                              partition_file = partition_file)
  
  # Write dataframe of information about constraint trees
  constraint_df_path <- paste0(dataset_constraint_tree_dir, prefix, "_constraint_tree_parameters.csv")
  write.csv(constraint_df, constraint_df_path, row.names = FALSE)
  
  # Return the constraint tree dataframe
  return(constraint_df)
} #end function





#### Estimating ML trees using a constraint tree ####
run.iqtree.with.constraint.tree <- function(output_prefix, alignment_path, constraint_tree_file,
                                            best_model = NA, free_rate_categories = NA, gamma = NA,
                                            iqtree_path = "iqtree2", num_bootstraps = NA, num_threads = 1, run.iqtree = FALSE,
                                            partitioned_check = FALSE, partition_file = NA){
  # Function estimate a constrained maximum likelihood tree: requires an alignment, a model, a constraint/guide tree, and the iqtree2 location
  
  
  ### Add partition file if present
  if (partitioned_check == FALSE){
    partition_call <- ""
  } else if (partitioned_check == TRUE){
    partition_call <- paste0(" -p ", partition_file, " ") # should this be spp? See doco: http://www.iqtree.org/doc/Command-Reference
  } 
  
  ### Set best_model as model for IQ-Tree run, including free-rate categories if necessary
  # Check whether the best model and free-rate categories are included
  best.model.provided = !is.na(best_model)
  rate.categories.provided = !is.na(free_rate_categories)
  # Assemble the call for the model
  if (best.model.provided == FALSE){
    # If no best_model specified for IQ-Tree, use model finder (-m MFP) command
    model_call = " -m MFP "
  } else if (best.model.provided == TRUE & rate.categories.provided == FALSE){
    # Best model provided, but no rate categories
    # Use only the best model in the IQ-Tree call
    model_call = paste0(" -m ", best_model, " ")
  } else if (best.model.provided == TRUE & rate.categories.provided == TRUE){
    # Both the best model and the rate categories are provided
    # Check the number of rate categories in the model and in the rate categories
    num_cats_rate_categories <- length(strsplit(free_rate_categories, ",")[[1]])/2
    num_cats_model <- as.numeric(gsub("R", "", grep("R", strsplit(best_model, "\\+")[[1]], value = T)))
    # Create a nice model with both the best model and the free rate category (weights and rates)
    model_call = paste0(" -m '", best_model, "{", free_rate_categories, "}'", " ")
  }
  
  ### Check for a gamma shape parameter (alpha) call
  # Determine whether the gamma variable is present (or is NA)
  if (is.na(gamma) == FALSE){
    # Make sure gamma is a character vector
    if (class(gamma) != "character"){
      gamma <- as.character(gamma)
    }
    # Check whether the gamma variable is an alpha parameter (by checking if there are any commas present)
    gamma_split <- strsplit(gamma, split = ",")[[1]]
    if (length(gamma) == length(gamma_split)){
      # There is only one gamma value: gamma here the gamma shape parameter
      # Strip any spaces from the gamma value
      gamma_clean <- gsub(" ", "", gamma)
      # Create an IQ-Tree command for gamma
      gamma_call <- paste0(" -a ",gamma_clean, " ")
    } else if (length(gamma) < length(gamma_split)){
      # There are multiple values within gamma: gamma here is a list of the rates, not an alpha parameter
      # Do not create an IQ-Tree command for gamma
      gamma_call <- ""
    }
  } else if (is.na(gamma) == TRUE){
    # Gamma variable is NA - do not create an IQ-Tree command for gamma
    gamma_call <- ""
  }
  
  ### Check for bootstraps call
  if (is.na(num_bootstraps) == TRUE){
    bootstraps_call <- ""
  } else if (is.na(num_bootstraps) == FALSE){
    bootstraps_call <- paste0(" -bb ", num_bootstraps, " ")
  }
  
  ### Add prefix (to label the hypothesis tree files) if present
  if (is.na(output_prefix) == TRUE){
    # If prefix is NA, do nothing
    prefix_call <- ""
  } else if (is.na(output_prefix) == FALSE){
    # If prefix is NA, add prefix to command line 
    prefix_call <- paste0(" -pre ", output_prefix, " ")
  } 
  
  ### Assemble the remaining parts of the command
  alignment_call <- paste0(" -s ", alignment_path, " ")
  constraint_tree_call <- paste0(" -g ", constraint_tree_file, " ")
  num_threads_call <- paste0(" -nt ", num_threads, " ")
  
  
  #### Collate iqtree command
  iqtree_call <- paste0(iqtree_path, alignment_call,  partition_call, model_call, gamma_call, bootstraps_call, constraint_tree_call, num_threads_call, prefix_call)
  # Remove any double or triple spaces from the call
  iqtree_call <- gsub("   ", " ", iqtree_call)
  iqtree_call <- gsub("  ", " ", iqtree_call)
  
  
  ### Run IQ-Tree, if desired
  if (run.iqtree == TRUE){
    # Run IQ-Tree
    system(iqtree_call)
  } # end if (run.iqtree == TRUE)
  
  ### Return the IQ-Tree command
  return(iqtree_call)
} # end function


run.one.constraint.tree <- function(index, constraint_df, run.iqtree = TRUE){
  # Quick function to take in a dataframe, take relevant variables, and call the run.iqtree.with.constraint.tree function
  
  # Identify row
  row <- constraint_df[index, ]
  
  # Feed row information into function call
  # Call function with lapply whether run.iqtree = TRUE or run.iqtree = FALSE:
  #   either way, want to run function to print iqtree command line
  iqtree_command <- run.iqtree.with.constraint.tree(output_prefix = row$hypothesis_tree_prefixes, alignment_path = row$alignment_path, constraint_tree_file = row$constraint_tree_paths,
                                                    best_model = row$best_model, free_rate_categories = row$estimated_rates, gamma = row$estimated_gamma,
                                                    iqtree_path = row$iqtree_path, num_bootstraps = row$num_bootstraps, num_threads = row$num_threads, run.iqtree = run.iqtree,
                                                    partitioned_check = FALSE, partition_file = row$partition_file)
  # Return the iqtree_command as output
  return(iqtree_command)
} # end function


run.one.constraint.dataframe <- function(csv_file, run.iqtree = TRUE){
  # Quick function to take in a dataframe, and estimate hypothesis trees by feeding it row by row into the run.one.constraint.tree function
  
  # Open the dataframe
  constraint_df <- read.csv(csv_file)
  # Estimate an ML tree in IQ-Tree for each constraint tree
  # Call function with lapply whether run.iqtree = TRUE or run.iqtree = FALSE:
  #   either way, want to run function to print iqtree command line
  lapply(1:nrow(constraint_df), run.one.constraint.tree, constraint_df = constraint_df, run.iqtree = run.iqtree)
} # end function





#### Collating multiple trees into a single file ####
combine.hypothesis.trees <- function(tree_id, constraint_tree_directory, outgroup_taxa = NA){
  # Function to open all hypothesis trees with a given id in a folder and collate them into one file
  
  # List all hypothesis trees in the folder
  all_constraint_tree_dir_files <- list.files(constraint_tree_directory, recursive = TRUE)
  # Remove any files with "ignore" in the name
  all_constraint_tree_dir_files <- grep("ignore", all_constraint_tree_dir_files, value = TRUE, invert = TRUE)
  # Find all files for this tree_id
  tree_id_files <- grep(tree_id, all_constraint_tree_dir_files, value = TRUE)
  # Find all hypothesis trees for this tree_id (hypothesis trees are marked by HX, where 1<= X <= 5)
  hypothesis_tree_files <- grep("H1|H2|H3|H4|H5", tree_id_files, value = TRUE)
  hypothesis_tree_treefiles <- grep("treefile", tree_id_files, value = TRUE)
  # Extend file path
  if (length(hypothesis_tree_treefiles) > 0){
    hypothesis_tree_treefiles <- paste0(constraint_tree_directory, hypothesis_tree_treefiles)
  }
  
  # Read in hypothesis tree files
  hypothesis_trees <- lapply(hypothesis_tree_treefiles, read.tree)
  # Convert hypothesis_trees from a list into a multiPhylo object 
  class(hypothesis_trees) <- "multiPhylo"
  
  # Output the (unrooted) hypothesis trees
  unrooted_file <- paste0(constraint_tree_directory, tree_id, "_unrooted_hypothesis_trees.tre")
  write.tree(hypothesis_trees, file = unrooted_file)
  
  if (class(outgroup_taxa) == "character"){
    # If the outgroup taxa are provided, root the hypothesis trees and save the rooted trees
    # Root hypothesis trees
    rooted_hypothesis_trees <- root(hypothesis_trees, outgroup_taxa)
    # Output the rooted hypothesis trees
    rooted_file <- paste0(constraint_tree_directory, tree_id, "_rooted_hypothesis_trees.tre")
    write.tree(hypothesis_trees, file = rooted_file)
    # Return paths for both rooted and unrooted hypothesis trees
    op_vec <- c(rooted_file, unrooted_file)
    names(op_vec) <- c("rooted_hypothesis_tree_file", "unrooted_hypothesis_tree_file")
  } else if (class(outgroup_taxa) == "logical"){
    # If the outgroup taxa are not provided, return only the path to the unrooted hypothesis trees
    op_vec <- c(unrooted_file)
    names(op_vec) <- c("unrooted_hypothesis_tree_file")
  }
  
  # Return the file names
  return(op_vec)
} # end function





#### Applying the mixture of trees model ####
tree.mixture.wrapper <- function(i, iqtree_tm_path, iqtree_num_threads, df){
  # Wrapper function to take the mixture of trees model 
  
  # Extract a row from the alignment
  row <- df[i, ]
  
  # Apply mixture of trees method with best model from maximum likelihood tree estimation
  # Run with +TR option (same branch lengths) 
  treemix_call <- run.tree.mixture.model(alignment_file = row$alignment_file, hypothesis_tree_file = row$hypothesis_tree_files, 
                                         partition_file = NA, use.partition = FALSE, prefix = paste0(row$prefix,".TR"),
                                         model = row$best_model, iqtree2_tree_mixtures_implementation = iqtree_tm_path, 
                                         tree_branch_option = "TR", number_parallel_threads = iqtree_num_threads,
                                         run.iqtree = FALSE)
  # Return the tree mixture iqtree call
  return(treemix_call)
}


run.tree.mixture.model <- function(alignment_file, hypothesis_tree_file, prefix, model, 
                                   iqtree2_tree_mixtures_implementation, tree_branch_option = "TR",
                                   number_parallel_threads, run.iqtree = TRUE){
  # Function runs the IQ-Tree2 mixture of trees model implementation given a sequence alignment, a set of hypothesis trees, and details about the model.
  # Currently cannot run with partition model
  # Tree branch options: 
  #   tree_branch_option = "TR" <- trees will have same length branches (for branches that appear in 2 or more trees)
  #   tree_branch_option = "T"  <- tree branches can have different lengths on different trees
  
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
    model_call = paste0("'", model_call, "+", tree_branch_option, "'")
  }
  
  # Add prefix if present
  if (is.na(prefix) == TRUE){
    # If prefix is NA, do nothing
    prefix_call <- ""
  } else if (is.na(prefix) == FALSE){
    # If prefix is NA, add prefix to command line 
    prefix_call <- paste0("-pre ", prefix)
    # Pad prefix_call with white space (for pasting into command line)
    prefix_call <- paste0(" ", prefix_call, " ")
  }
  
  # Assemble the command for the tree mixtures model
  treemix_command <- paste0(iqtree2_tree_mixtures_implementation, " -s ", alignment_file, 
                            " -m  ", model_call, " -te ", hypothesis_tree_file, 
                            " -nt ", number_parallel_threads, prefix_call)
  
  # Change working directories (to store IQ-Tree output files in the right place)
  setwd(dirname(hypothesis_tree_file))
  
  if (run.iqtree == TRUE){
    # Call IQ-Tree2 with the command
    system(treemix_command)
  } # end if (run.iqtree == TRUE)
  
  
  # Return the iqtree command line
  return(treemix_command)
} # end function





