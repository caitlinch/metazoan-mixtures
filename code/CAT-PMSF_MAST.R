### caitlinch/metazoan-mixtures/code/CAT-PMSF_MAST.R
# This script performs analyses using the CAT-PMSF model for 14 empirical data sets
# Caitlin Cherryh 2025

# This script is designed to be used in conjunction with the Rproject for this repository

## TODO
# - Fixed topology used to infer site-specific amino acid frequency profiles -
#   should this be a bifucating tree topology with branch lengths? If so, I will
#   need to estimate constrained trees for each dataset under the CAT-POISSON model.
# - Need to repeat process for each topology under consideration - i.e.,
#   CAT-POISSON-PORI, CAT-POISSON-CTEN. This doubles the number of necessary
#   analyses and for each dataset (at minimum) I would need to run:
#           - MAST 2-tree CAT-POISSON-CTEN
#           - MAST 2-tree CAT-POISSON-PORI
#           - MAST 5-tree CAT-POISSON-CTEN
#           - MAST 5-tree CAT-POISSON-PORI
# - I wasn't able to install phylobayes-mpi on Dayhoff (issue with MPI compilers)



#### 01. Load packages ####
library(yaml)

source("code/func_CAT-PMSF.R")



#### 02. Input parameters ####
# Open config file to get file paths
config <- read_yaml("code/CAT-PMSF_config.yaml")
params <- config[[which(names(config) == config$environment)]]
# Open the dataframe with the alignment details
alignment_df <- read.csv("output/alignment_dimensions.csv")
alignment_df$full_path <- paste0(params$alignments, list.files(params$alignments))
# Assemble constraint tree list
constraint_trees <- paste0(params$constraint_trees, list.files(params$constraint_trees))
# Filter constraint tree list to include only desired topologies
profile_topologies <- grep("constraint_tree_1|constraint_tree_2", constraint_trees, value = T)



#### 03. PhyloBayes model inference ####
# Following the procedure from Giacomelli et al. 2025 (Genome Biol. Evol.
#     17:evae273, doi:10.1093/gbe/evae273) and the associated GitHub repository
#     (https://github.com/mgiacom/tardigrades_catpmsf/)
# Based on command lines from the file:
#     https://github.com/mgiacom/tardigrades_catpmsf/blob/main/pipeline
# Steps:
#   - Infer compositional profiles under CAT-POISSON in PhyloBayes
#   - Check parameter convergence
#   - Estimate site-specific profiles under CAT-POISSON
#   - Convert PhyloBayes .siteprofiles files to IQ-Tree .sitefreq file format
#         using script from convert-site-dists.py from Szánthó et al. 2023
#         (Sys. Biol. 72:767–780, doi:10.1093/sysbio/syad013) and the associated
#         GitHub repo (https://github.com/drenal/cat-pmsf-paper)

## Create command lines to infer PhyloBayes command lines
# pb_mpi -s matrix_40sp.fasta -T sister_to_lobopodia.tree -cat -poisson sister_to_lobopodia_chain1 (and chain_2)
# tracecomp -x 500 sister_to_lobopodia_chain1 sister_to_lobopodia_chain2
# readpb_mpi -ss -x 500 10 sister_to_lobopodia_chain1

# Convert to IQ-Tree file format: convert-site-dists.py sister_to_lobopodia_chain1.siteprofiles

#### 04. Unconstrained tree estimation ####




#### 05. Constrained tree estimation ####




#### 06. MAST model estimation ####




#### 07. AU test ####



