## caitlinch/metazoan-mixtures/code/func_BIC.R
# Functions for manipulating, processing, and preparing datasets
# Caitlin Cherryh 2023


#### Packages ####
library(ape) # read.tree
library(phangorn) # as.splits
library(dplyr) # manipulating dataframes


#### Calculating number of different splits ####
compare.splits.2.trees <- function(trees_path){
  ## Take 2 trees and calculate the number of different splits
  # Calculate the number of different splits using the RF distance
  t1 <- read.tree(trees_path[1])
  t2 <- read.tree(trees_path[2])
  num_different_splits <- RF.dist(t1, t2)
  return(num_different_splits)
}




