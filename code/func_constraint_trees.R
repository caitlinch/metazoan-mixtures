## caitlinch/metazoan-mixtures/func_constraint_trees.R
# Caitlin Cherryh 2022

# Functions for testing and applying constraint trees in iqtree2

run.iqtree.with.constraint.tree <- function(alignment_path, constraint_tree_file, iqtree_path, model = NA, num_threads = 1, prefix = NA){
  # Function to apply IQ-Tree to a series of alignments with a constraint tree
  
  # Set model for IQ-Tree run
  if (is.na(model) == TRUE){
    # If no model specified for IQ-Tree, use model finder (-m MFP) command
    iq_model = "MFP"
  } else {
    # Otherwise, use model specified in function call
    iq_model = model
  }
  
  # Create command line for IQ-Tree2
  if (is.na(prefix) == TRUE){
    iqtree_call <- paste0(iqtree_path, " -s ", alignment_path, " -m ", iq_model, " -g ", constraint_tree_file, " -T ", num_threads)
  } else (is.na(prefix) == FALSE){
    iqtree_call <- paste0(iqtree_path, " -s ", alignment_path, " -m ", iq_model, " -g ", constraint_tree_file, " -T ", num_threads,  " --prefix ", prefix)
  }
  
  # Run IQ-Tree
  system(iqtree_call)
}



apply.one.constraint.tree <- function(index, df){
  # Quick function to take in a dataframe, take relevant variables, and call the run.iqtree.with.constraint.tree function
  
  # Identify row
  row <- df[index, ]
  
  # Feed row information into function call
  run.iqtree.with.constraint.tree(alignment_path = row$alignment_path, constraint_tree_file = row$constraint_tree_paths, 
                                  iqtree_path = row$iqtree_path, prefix = row$constraint_prefixes, model = row$model)
  
}