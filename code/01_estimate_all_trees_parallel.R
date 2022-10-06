# metazoan-mixtures/code/01_estimate_all_trees_parallel.R
## This script estimates maximum likelihood trees and maximum likelihood under constraint trees for empirical data sets
# Caitlin Cherryh, 2022

## This script:
# 1. Estimates ML trees for empirical data sets under different models
# 2. Estimates ML trees using a constraint tree for empirical data sets under different models
# 3. Applies the MAST (Mixtures Across Sites and Trees) model 



#### 1. Input parameters ####
## Specify parameters:
# alignment_dir       <- Directory containing alignments for all data sets
#                        Alignments have the naming convention dataset.matrix_name.sequence_type.fa
#                        E.g. Cherryh2022.all_taxa.aa.fa
# output_dir          <- Directory for IQ-Tree output (trees and tree mixtures)
# repo_dir            <- Location of caitlinch/metazoan-mixtures github repository

# iqtree2             <- Location of IQ-Tree2 stable release
# iqtree_tm           <- Location of IQ-Tree2 MAST release

# iqtree_num_threads  <- Number of parallel threads for IQ-Tree to use. Can be a set number (e.g. 2) or "AUTO"

location = "dayhoff"
if (location == "local"){
  alignment_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/01_Data_all/"
  output_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/"
  repo_dir <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"
  
  iqtree2 <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_Software_IQ-Tree/iqtree-2.2.0-MacOSX/bin/iqtree2"
  iqtree2_tm <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_Software_IQ-Tree/iqtree-2.2.0.7.mix-MacOSX/bin/iqtree2"
  
  iqtree_mrate <- "E,I,G,I+G,R,I+R"
  iqtree_num_threads <- 1
  ml_tree_bootstraps <- 1000
  parallel_threads <- 1
  
} else if (location == "soma"){
  alignment_dir <- "/data/caitlin/metazoan-mixtures/data_all/"
  output_dir <- "/data/caitlin/metazoan-mixtures/output/"
  repo_dir <- "/data/caitlin/metazoan-mixtures/"
  
  iqtree2 <- "/data/caitlin/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2"
  iqtree2_tm <- "/data/caitlin/metazoan-mixtures/iqtree/iqtree-2.2.0.7.mix-Linux/bin/iqtree2"
  
  iqtree_mrate <- "E,I,G,I+G,R,I+R"
  iqtree_num_threads <- 1
  ml_tree_bootstraps <- 1000
  parallel_threads <- 20
  
} else if (location == "dayhoff"){
  alignment_dir <- "/home/u5348329/metazoan-mixtures/data_all/"
  output_dir <- "/home/u5348329/metazoan-mixtures/output/"
  repo_dir <- "/home/u5348329/metazoan-mixtures/"
  
  iqtree2 <- "/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2"
  iqtree2_tm <- "/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0.7.mix-Linux/bin/iqtree2"
  
  iqtree_mrate <- "E,I,G,I+G,R,I+R"
  iqtree_num_threads <- 1
  ml_tree_bootstraps <- 1000
  parallel_threads <- 20
} else if (location == "laptop"){
  alignment_dir <- "/Users/caitlin/Documents/PhD/Ch03_sponge_mixtures/01_alignments/"
  output_dir <- "/Users/caitlin/Documents/PhD/Ch03_sponge_mixtures/02_output/"
  repo_dir <- "/Users/caitlin/Repositories/metazoan-mixtures/"
  
  iqtree2 <- "/Users/caitlin/Documents/PhD/Ch03_sponge_mixtures/iqtree-2.2.0-MacOSX/bin/iqtree2"
  iqtree2_tm <- "/Users/caitlin/Documents/PhD/Ch03_sponge_mixtures/iqtree-2.2.0.7.mix-MacOSX/bin/iqtree2"
  
  iqtree_mrate <- "E,I,G,I+G,R,I+R"
  iqtree_num_threads <- 1
  ml_tree_bootstraps <- 1000
  parallel_threads <- 1
}



#### 2. Prepare functions, variables and packages ####
# Open packages
library(parallel)

# Source functions
source(paste0(repo_dir, "code/func_constraint_trees.R"))
source(paste0(repo_dir, "code/func_data_processing.R"))

# Source information about datasets
source(paste0(repo_dir, "code/data_dataset_info.R"))

# Extract the list of all files from the folder containing alignments/models/partition files
all_files <- list.files(alignment_dir)
if (length(all_files) > 0){
  all_files <- paste0(alignment_dir, all_files)
}
# Extract the list of alignments (i.e. files that contain the word "alignment")
all_alignments <- grep("\\.alignment\\.", all_files, value = T)

# Sort and categorise models of sequence evolution 
all_models <- sort(unique(unlist(lapply(all_models, sort.model.chunks))))
# Identify unique matrices within the models
model_components <- sort(unique(unlist(strsplit(all_models, "\\+"))))
# Remove components that do not correspond to a matrix or model (e.g. remove "I", "F", "G4")
model_components <- model_components[!model_components %in% c("F", "FO", "G", "G4", "I", "R", "R4")]
# PMSF is more complex than other models - requires separate function
# For now, remove PMSF from list
model_components <- model_components[!model_components == "PMSF"]
# Move ModelFinder to the end of the list
model_components <- model_components[!model_components == "ModelFinder"]
model_components <- c(model_components, "ModelFinder")
# Note: partitioning schemes currently not possible in mixture of trees implementation



#### 3. Estimate maximum likelihood trees for each combination of model and dataset ####
# Create a folder for the ml trees and move to that folder
ml_tree_dir <- paste0(output_dir, "maximum_likelihood_trees/")
if (file.exists(ml_tree_dir) == FALSE){dir.create(ml_tree_dir)}
setwd(ml_tree_dir)

# Create a dataframe of combinations of alignments and models
ml_tree_df <- expand.grid(dataset = unlist(lapply(strsplit(basename(all_alignments), "\\."), "[[", 1)),
                          model_code = model_components,
                          stringsAsFactors = FALSE)
ml_tree_df$model_mset <- ml_tree_df$model_code
ml_tree_df$model_m <- NA
ml_tree_df$model_mrate <- iqtree_mrate
ml_tree_df$sequence_format = unlist(lapply(strsplit(basename(all_alignments), "\\."), "[[", 3))
ml_tree_df$matrix_name <- unlist(lapply(strsplit(basename(all_alignments), "\\."), "[[", 2))
ml_tree_df$prefix <- paste0(ml_tree_df$dataset, ".", ml_tree_df$matrix_name, ".", ml_tree_df$model_code)
ml_tree_df$iqtree_num_threads <- iqtree_num_threads
ml_tree_df$iqtree_num_bootstraps <- ml_tree_bootstraps
ml_tree_df$alignment_file <- all_alignments
ml_tree_df$iqtree_file <- paste0(ml_tree_df$prefix, ".iqtree")
ml_tree_df$ml_tree_file <- paste0(ml_tree_df$prefix, ".treefile")

# Fix model specification for rows with ModelFinder
ml_tree_df[ml_tree_df$model_mset == "ModelFinder",]$model_m <- "MFP"
ml_tree_df$model_mset[which(ml_tree_df$model_code == "ModelFinder")] <- NA

# Sort matrix by dataset and matrix
ml_tree_df <- ml_tree_df[order(ml_tree_df$dataset, ml_tree_df$matrix_name),]

# Create IQ-Tree2 commands
ml_tree_df$iqtree2_call <- unlist(lapply(1:nrow(ml_tree_df), ml.iqtree.wrapper, iqtree_path = iqtree2, df = ml_tree_df))

# Save dataframe
ml_tree_df_name <- paste0(output_dir, "maximum_likelihood_tree_estimation_parameters.csv")
write.csv(ml_tree_df, file = ml_tree_df_name, row.names = FALSE)

# Run IQ-Tree commands to estimate ML trees for each model/matrix combination
mclapply(ml_tree_df$iqtree2_call, system, mc.cores = parallel_threads)



#### 4. Estimate constraint and hypothesis trees for each combination of model and dataset ####
# Create a folder for the ml trees and move to that folder
c_tree_dir <- paste0(output_dir, "constraint_trees/")
if (file.exists(c_tree_dir) == FALSE){dir.create(c_tree_dir)}
setwd(c_tree_dir)

# Extract the best model for each combination of matrix and model
ml_tree_df$best_model <- unlist(lapply(ml_tree_df$iqtree_file, extract.best.model))
# Save dataframe
ml_tree_df_name <- paste0(output_dir, "maximum_likelihood_tree_estimation_parameters_complete.csv")
write.csv(ml_tree_df, file = ml_tree_df_name, row.names = FALSE)

# Create a constraint df for each row in the ml_tree_df
constraint_list <- lapply(1:nrow(ml_tree_df), constraint.tree.wrapper, output_directory = c_tree_dir, 
                          iqtree_path = iqtree2, iqtree_num_threads = iqtree_num_threads, 
                          dataset_info = all_datasets, matrix_taxa_info = matrix_taxa,
                          df = ml_tree_df)
# Collate the constraints into a single dataframe
constraint_df <- do.call(rbind, constraint_list)
# Add the mrate = NA options for IQ-Tree to the dataframe (do not include mrate option for estimating constraint trees)
constraint_df$model_mrate <- NA

# Estimate hypothesis trees for each of the constraint trees (call one row of the dataframe at a time)
constraint_df$iqtree2_call <- unlist(lapply(1:nrow(constraint_df), run.one.constraint.tree, df = constraint_df, run.iqtree = FALSE))

# Save the constraint tree dataframe
c_tree_df_name <- paste0(output_dir, "constraint_tree_estimation_parameters.csv")
write.csv(constraint_df, file = c_tree_df_name, row.names = FALSE)

# Run IQ-Tree commands to estimate hypothesis trees for each model/matrix combination
mclapply(ml_tree_df$iqtree2_call, system, mc.cores = parallel_threads)

# Combine hypothesis trees into one file and save
ml_tree_df$hypothesis_tree_files <- lapply(ml_tree_df$prefix, combine.hypothesis.trees, constraint_tree_directory = c_tree_dir, 
                                           outgroup_taxa = NA)
# Save dataframe
ml_tree_df_name <- paste0(output_dir, "MAST_estimation_parameters.csv")
write.csv(ml_tree_df, file = ml_tree_df_name, row.names = FALSE)



#### 5. Apply mixtures across trees and sites (MAST model) ####
# Create a folder for the ml trees and move to that folder
m_tree_dir <- paste0(output_dir, "tree_mixtures/")
if (file.exists(m_tree_dir) == FALSE){dir.create(m_tree_dir)}
setwd(m_tree_dir)

# Create the tree mixture prefix and command lines
ml_tree_df$MAST_prefix <- paste0(ml_tree_df$prefix,".TR")
ml_tree_df$MAST_call <- lapply(1:nrow(ml_tree_df), tree.mixture.wrapper, iqtree_tm_path = iqtree2_tm, 
                               iqtree_num_threads = iqtree_num_threads, df = ml_tree_df)

# Run the mixture of trees models
mclapply(ml_tree_df$MAST_call, system, mc.cores = parallel_threads)


############### Incomplete code: still need to extract results from MAST model and save to csv file
# # Identify iqtree files from tree mixture run
# all_files <- paste0(a_tm_op_dir, list.files(a_tm_op_dir))
# all_iqtree_files <- grep("\\.iqtree", all_files, value = TRUE)
# tree_mixture_tr_iqfile <- grep(paste0(a_m_prefix,".TR"), all_iqtree_files, value = TRUE)
# 
# # Extract information about each tree mixture model run
# tr_results <- extract.tree.mixture.results(tree_mixture_file = tree_mixture_tr_iqfile, 
#                                            dataset = a_dataset, prefix = paste0(a_m_prefix,".TR"), 
#                                            model = m, best_model = best_model, tree_branch_option = "TR")
# 
# # Output results dataframe
# op_file <- paste0(a_tm_op_dir, a_m_prefix, "_tree_mixture_results.csv")
# write.csv(tr_results, file = op_file, row.names = FALSE)

