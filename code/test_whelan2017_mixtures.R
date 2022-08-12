## caitlinch/metazoan-mixtures/test_whelan2017_mixtures.R
# Caitlin Cherryh 2022

# Proof of concept for the mixture of trees method

#### Step 1: Input parameters ####
# main_dir                <- path to caitlinch/metazoan-mixtures git repository
# gene_folder             <- path to folder containing fasta files for each gene in the Whelan 2017 dataset
# iqtree_path             <- path to IQ-Tree2 executable with mixtures of trees implementation
# constraint_tree_dir     <- folder to store constraint trees in
# number_parallel_threads <- number of cores to use for parallel processes

location = "local"
if (location == "local"){
  main_dir <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"
  data_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/02_Data_processed/"
  constraint_tree_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_hypothesis_trees/"
  
  number_parallel_threads = "AUTO"
} else if (location == "soma"){
  main_dir <- "/data/caitlin/metazoan-mixtures/"
  gene_folder <- "/data/caitlin/metazoan-mixtures/data_whelan2017/genes/"
  iqtree_path <- "data/caitlin/metazoan-mixtures/iqtree/iqtree-2.2.0.6.mix-Linux/bin/iqtree2"
  constraint_tree_dir <- "/data/caitlin/metazoan-mixtures/constraint_trees/"
  
  number_parallel_threads = 20
}

assemble_constraint_trees <- FALSE
estimate_constraint_trees <- FALSE
apply_tree_mixtures <- TRUE



#### Step 2: Prepare analysis ####
# Source function files
source(paste0(main_dir, "code/func_constraint_trees.R"))

# Open packages
library(parallel)
library(ape)
library(phangorn)

# Create folders if necessary
if (dir.exists(constraint_tree_dir) == FALSE){dir.create(constraint_tree_dir)}



#### Step 3: Prepare constraint trees ####
if (assemble_constraint_trees == TRUE){
  ## For Whelan2017 data:
  # Set dataset name
  dataset = "Whelan2017"
  # Select model of sequence evolution 
  model = "CAT"
  
  # Create folder for each dataset inside the constraint tree folder
  dataset_constraint_tree_dir <- paste0(constraint_tree_dir, dataset, "/")
  if (dir.exists(dataset_constraint_tree_dir) == FALSE){dir.create(dataset_constraint_tree_dir)}
  
  
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
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, dataset, "_constraint_tree_", "1", ".nex")
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
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, dataset, "_constraint_tree_", "2", ".nex")
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
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, dataset, "_constraint_tree_", "3", ".nex")
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
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, dataset, "_constraint_tree_", "4", ".nex")
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
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, dataset, "_constraint_tree_", "5", ".nex")
  write(constraint_tree_5, file = constraint_tree_file_name)
  
  # Assemble dataframe of information about the constraint trees
  constraint_df <- data.frame(constraint_tree_id = 1:5,
                              constraint_tree_paths = paste0(dataset_constraint_tree_dir, dataset, "_constraint_tree_", 1:5, ".nex"),
                              constraint_prefixes = paste0(dataset, "_", model, "_ML_H", 1:5),
                              alignment_path = gene_folder,
                              model = "CAT",
                              iqtree_path = iqtree_path,
                              constraint_trees = c(constraint_tree_1, constraint_tree_2, constraint_tree_3, 
                                                   constraint_tree_4, constraint_tree_5),
                              num_threads = number_parallel_threads)
  
  # Write dataset of information about constraint trees
  constraint_df_path <- paste0(dataset_constraint_tree_dir, dataset, "_constraint_tree_parameters.csv")
  write.csv(constraint_df, constraint_df_path, row.names = FALSE)
}



#### Step 4: Estimate trees with constraint trees ####
if (estimate_constraint_trees == TRUE){
  ## For Whelan2017 data:
  # Set dataset name
  dataset = "Whelan2017"
  
  # Create folder for each dataset inside the constraint tree folder
  dataset_constraint_tree_dir <- paste0(constraint_tree_dir, dataset, "/")
  if (dir.exists(dataset_constraint_tree_dir) == FALSE){dir.create(dataset_constraint_tree_dir)}
  
  
  # Set working directory to dataset_constraint_tree_dir so IQ-Tree output is saved with the constraint trees
  setwd(dataset_constraint_tree_dir)
  
  # For trees with all 76 taxa
  estimate_trees_df <- read.csv(paste0(dataset_constraint_tree_dir, dataset, "_constraint_tree_parameters.csv"))
  # Estimate an ML tree in IQ-Tree for each constraint tree
  lapply(1:nrow(estimate_trees_df), apply.one.constraint.tree, estimate_trees_df)
}



#### Step 5: Collate trees ####
# Set dataset
dataset <- "Whelan2017"
# List all files in the constraint tree directory
all_dir_files <- list.files(constraint_tree_dir, recursive = TRUE)
# Remove previous attempts ("Old") or test runs ("Small")
all_dir_files <- grep("Old", all_dir_files, value = TRUE, invert = TRUE)
all_dir_files <- grep("Small", all_dir_files, value = TRUE, invert = TRUE)
# Identify constraint tree files for the specified dataset
dataset_files <- grep(dataset, all_dir_files, value = TRUE)
constraint_tree_files <- grep(".treefile", dataset_files, value = TRUE)
constraint_trees <- grep(".treefile", constraint_tree_files, value = TRUE)
# Open constraint trees
ctrees <- read.tree(text = unlist(lapply(paste0(constraint_tree_dir, constraint_trees), readLines))) 
rooted_ctrees <- root(ctrees, c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"))
# Write constraint trees
treemix_tree_file <- paste0(constraint_tree_dir, dataset, "_hypothesis_trees.tre")
write.tree(rooted_ctrees, file = treemix_tree_file)



#### Step 6: Apply the mixture of trees method ####
if (apply_tree_mixtures == TRUE){
  treemix_command <- paste0(iqtree_path, " -s ", gene_folder, " -m  'LG+TR' -te ", treemix_tree_file, " -nt 20 -pre Whelan2017_LG_5TR_1Q")
  system(treemix_command)
}



