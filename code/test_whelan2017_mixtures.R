## caitlinch/metazoan-mixtures/test_whelan2017_mixtures.R
# Caitlin Cherryh 2022

# Proof of concept for the mixture of trees method

#### Step 1: Input parameters ####
# main_dir                <- path to caitlinch/metazoan-mixtures git repository
# data_dir                <- path to folder containing concatenated alignments for all datasets  (and model/partition files if necessary)
# constraint_tree_dir     <- folder to store constraint trees in
# iqtree2                 <- path to IQ-Tree2 stable implementation (v2.2.0 COVID-edition)
# iqtree2_tm              <- path to IQ-Tree2 executable with mixtures of trees implementation (v2.2.0.6.mix)
# number_parallel_threads <- number of cores to use for parallel processes

# datasets                <- List of identifiers for datasets (first author surname + year of publication)
# datasets_to_run         <- List of dataset identifiers to run through pipeline. Can contain all datasets (if datasets_to_run = datasets), or a selected subset

# assemble_constraint_trees  <- Whether to assemble the constraint trees using the taxa from each dataset (TRUE or FALSE)
# estimate_hypothesis_trees  <- Whether to run IQ-Tree with the constraint trees to estimate an ML tree for each constraint (TRUE or FALSE)
# apply_tree_mixtures        <- Whether to apply the mixture of trees model to the hypothesis trees (TRUE or FALSE)

location = "local"
if (location == "local"){
  main_dir <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"
  data_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/02_Data_processed/"
  constraint_tree_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_hypothesis_trees/"
  iqtree2 <- "iqtree2"
  iqtree2_tm <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_Software_IQ-Tree/iqtree-2.2.0.6.mix/iqtree-2.2.0.6.mix-MacOSX/bin/iqtree2"
  number_parallel_threads = "AUTO"
} else if (location == "soma"){
  main_dir <- "/data/caitlin/metazoan-mixtures/"
  data_dir <- "/data/caitlin/metazoan-mixtures/data_all/"
  constraint_tree_dir <- "/data/caitlin/metazoan-mixtures/constraint_trees/"
  iqtree2 <- "/data/caitlin/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2"
  iqtree2_tm <- "data/caitlin/metazoan-mixtures/iqtree/iqtree-2.2.0.6.mix-Linux/bin/iqtree2"
  number_parallel_threads = 20
}

datasets <- c("Whelan2017")
datasets_to_run <- datasets

assemble_constraint_trees <- FALSE
estimate_hypothesis_trees <- FALSE
apply_tree_mixtures <- TRUE

# Test parameters
dataset = "Whelan2017"
model = "LG"



#### Step 2: Prepare analysis ####
# Source function files
source(paste0(main_dir, "code/func_constraint_trees.R"))

# Source clade and taxa file
source(paste0(main_dir, "code/data_dataset_info.R"))

# Open packages
library(parallel)
library(ape)
library(phangorn)

# Create folders if necessary
if (dir.exists(constraint_tree_dir) == FALSE){dir.create(constraint_tree_dir)}



#### Step 3: Prepare constraint trees ####
if (assemble_constraint_trees == TRUE){
  for (dataset in datasets_to_run){
    
    # Create folder for each dataset inside the constraint tree folder
    dataset_constraint_tree_dir <- paste0(constraint_tree_dir, dataset, "/")
    if (dir.exists(dataset_constraint_tree_dir) == FALSE){dir.create(dataset_constraint_tree_dir)}
    
    #Extract the list of information about this dataset
    dataset_list <- all_datasets[[dataset]]
    # Separate out the relevant groups of taxa from the dataset of interest
    outgroup_taxa = dataset_list$Outgroup
    ctenophora_taxa = dataset_list$Ctenophora
    porifera_taxa = dataset_list$Porifera
    placozoa_taxa = dataset_list$Placozoa
    cnidaria_taxa = dataset_list$Cnidaria
    bilateria_taxa = dataset_list$Bilateria
    # Identify the two groups of sponges
    sponges_1_taxa <- unlist(dataset_list[c(dataset_list$Sponges_1)])
    sponges_2_taxa <- unlist(dataset_list[c(dataset_list$Sponges_2)])
    # Identify model
    model <- dataset_list$Models
    
    # Set ModelFinder to run if ModelFinder or PartitionFinder were used in the original study
    if (model == "PartitionFinder" | model == "ModelFinder"){
      model = "MFP+MERGE -spp"
      model_id = "ModelFinder"
    }
    
    # Create the constraint trees and constraint tree information dataframe
    constraint_df <- create.constraint.trees(dataset, dataset_constraint_tree_dir, model, model_id, outgroup_taxa, ctenophora_taxa, porifera_taxa,
                                             sponges_1_taxa, sponges_2_taxa, placozoa_taxa, cnidaria_taxa, bilateria_taxa)
  } # end for (dataset in datasets_to_run)
} # end if (assemble_constraint_trees == TRUE)



#### Step 4: Estimate trees with constraint trees ####
if (estimate_hypothesis_trees == TRUE){
  for (dataset in datasets_to_run){
    
    # Create folder for each dataset inside the constraint tree folder (if it doesn't already exist)
    dataset_constraint_tree_dir <- paste0(constraint_tree_dir, dataset, "/")
    if (dir.exists(dataset_constraint_tree_dir) == FALSE){dir.create(dataset_constraint_tree_dir)}
    
    # Set working directory to dataset_constraint_tree_dir so IQ-Tree output is saved with the constraint trees
    setwd(dataset_constraint_tree_dir)
    
    # Find a list of the constraint dataframe files
    all_param_files <- grep(".csv",list.files(dataset_constraint_tree_dir), value = TRUE)
    # Select the constraint tree files for this dataset
    dataset_param_files <- grep(dataset, all_param_files, value = TRUE)
    # Run one dataset at a time
    if (length(dataset_param_files) > 0){
      lapply(paste0(dataset_constraint_tree_dir, dataset_param_files), run.one.constraint.dataframe)
    } # end if (length(dataset_param_files) > 0)
    
  } # end for (dataset in datasets_to_run)
} # end if (estimate_hypothesis_trees == TRUE)



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



