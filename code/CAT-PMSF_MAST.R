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
## Set file paths
# Open config file to get file paths
params <- read_yaml("code/CAT-PMSF_config.yaml")
params <- params[[which(names(params) == params$environment)]]
# Add output paths for specific analyses
params$output_profiles <- paste0(params$output, "01_CAT-PMSF_profiles/")
params$output_trees <- paste0(params$output, "02_unconstrained_trees/")
params$output_constrained_trees <- paste0(params$output, "03_constrained_trees/")
params$output_MAST <- paste0(params$output, "04_MAST/")
params$output_AU_test <- paste0(params$output, "04_AU_test/")

## Alignment details
# Open the data frame with the alignment details
alignment_df <- read.csv("output/alignment_dimensions.csv")
alignment_df$full_path <- paste0(params$alignments, alignment_df$alignment_path)

## Constraint trees
# Filter constraint tree list to include only desired topologies
#   Models included in constraint trees are: LG4M, UL3, PMSF LG+C60, PMSF C60, LG+C20, LG+C60, C60
# Want only constrained LG+C20, LG+C60 or C60 models for CTEN-sister (ML_H1) and PORI-sister (ML_H2)
profile_topologies <-
  grep("LG_C20|LG_C60|C60",
       grep("ML_H1|ML_H2",
            grep("\\.treefile",
                 grep("PMSF",
                      list.files(params$constrained_CXX_trees),
                      invert = T, value = T),
                 value = T),
            value = T),
       value = T)



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
# Remember: these commands should be run in the params$output_profiles directory
commands_df <- as.data.frame(
  do.call(rbind,
          lapply(1:nrow(alignment_df),
                 infer.phylobayes.profile.wrapper,
                 alignment_df, params, profile_topologies)
  )
)

## Create command lines to convert site profiles to IQ-Tree file format (.sitefreq)
commands_df <- as.data.frame(
  do.call(rbind,
          lapply(
            1:nrow(commands_df),
            convert.profile.for.iqtree,
            commands_df, params)
  )
)

## Write csv file
write.csv(
  commands_df,
  file = paste0(params$output, "CAT-PMSF_phylobayes_command_lines.csv")
)

#### 04. Unconstrained tree estimation ####




#### 05. Constrained tree estimation ####




#### 06. MAST model estimation ####




#### 07. AU test ####



