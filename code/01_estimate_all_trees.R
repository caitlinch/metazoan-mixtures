# metazoan-mixtures/code/01_estimate_all_trees.R
## This script estimates maximum likelihood trees and maximum likelihood under constraint trees for empirical datasets
# Caitlin Cherryh, 2022

## This script:
# 1. Estimates ML trees for empirical datasets under different models
# 2. Estimates ML trees using a constraint tree for empirical datasets under different models
# 3. Applies the MAST (Mixtures Across Sites and Trees) model 



#### 1. Parameters ####
## Specify directories:



# - Take a single alignment:
#   - Take a single model of sequence evolution
#     - Use ModelFinder to determine the best model within that model category
#     - Estimate an ML tree with that model of sequuence evolution
#     - Estimate the constraint trees
#     - Apply the mixture of trees method using +TR
#     - Apply the mixture of trees method using +T