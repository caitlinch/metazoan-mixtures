## caitlinch/metazoan-mixtures/func_constraint_trees.R
# Caitlin Cherryh 2022

# Functions for testing and applying constraint trees in iqtree2

run.iqtree.with.constraint.tree <- function(alignment_path, constraint_tree_file, iqtree_path, prefix, model = NA){
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
  iqtree_call <- paste0(iqtree_path, " -s ", alignment_path, " -m ", iq_model, " -g ", constraint_tree_file, " --prefix ", prefix)
  
  # Run IQ-Tree
  system(iqtree_call)
}