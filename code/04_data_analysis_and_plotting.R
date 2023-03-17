# metazoan-mixtures/code/03_data_analysis_and_plotting.R
## This script performs data analysis and creates plots of methods and results
# Caitlin Cherryh, 2022

## This script:
# 1. Plots figures for the introduction and methods section of the manuscript



#### 1. Input parameters ####
## Specify parameters:
# output_dir          <- Directory for IQ-Tree output (trees and tree mixtures)
# results_dir         <- Directory for results and plots
# repo_dir            <- Location of caitlinch/metazoan-mixtures github repository

output_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/05_plotting/"
repo_dir <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"



#### 2. Prepare functions, variables and packages ####
# Open packages
library(ggplot2)
library(ggtree)
library(patchwork)

# Open function files
source(paste0(repo_dir,"code/func_plotting.R"))



#### 3. Plot figures for introduction and methods sections of manuscript ####
## Plotting the alternative phylogenetic hypotheses
# Open the trees
trees_file <- paste0(repo_dir, "trees/alternative_phylogenetic_hypotheses.nex")
trees <- read.tree(trees_file)

