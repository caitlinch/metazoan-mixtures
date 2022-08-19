# metazoan-mixtures/code/01_estimate_all_trees.R
## This script estimates maximum likelihood trees and maximum likelihood under constraint trees for empirical datasets
# Caitlin Cherryh, 2022

## This script:
# 1. Estimates ML trees for empirical datasets under different models
# 2. Estimates ML trees using a constraint tree for empirical datasets under different models
# 3. Applies the MAST (Mixtures Across Sites and Trees) model 



#### 1. Parameters ####
## Specify parameters:
# alignment_dir     <- Directory containing alignments for all datasets
# output_dir        <- Directory for IQ-Tree output (trees and tree mixtures)
# repo_dir          <- Location of caitlinch/metazoan-mixtures github repository

# iqtree2           <- Location of IQ-Tree2 stable release
# iqtree_tm         <- Location of IQ-Tree2 MAST release

location = "soma"
if (location == "local"){
  alignment_dir <- "/Users/caitlin/Documents/PhD/Ch03_sponge_mixtures/01_alignments/"
  output_dir <- "/Users/caitlin/Documents/PhD/Ch03_sponge_mixtures/02_output/"
  repo_dir <- "/Users/caitlin/Repositories/metazoan-mixtures/"
  
  iqtree2 <- "/Users/caitlin/Documents/PhD/Ch03_sponge_mixtures/iqtree-2.2.0-MacOSX/bin/iqtree2"
  iqtree2_tm <- "/Users/caitlin/Documents/PhD/Ch03_sponge_mixtures/iqtree-2.2.0.7.mix-MacOSX/bin/iqtree2"
  
} if (location == "soma"){
  alignment_dir <- ""
  output_dir <- ""
  repo_dir <- ""
}


#### 2. Process each dataset for each model of sequence evolution ####
# - Take a single alignment:
#   - Take a single model of sequence evolution
#     - Use ModelFinder to determine the best model within that model category
#     - Estimate an ML tree with that model of sequuence evolution
#     - Estimate the constraint trees
#     - Apply the mixture of trees method using +TR
#     - Apply the mixture of trees method using +T