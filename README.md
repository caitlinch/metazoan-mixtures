# metazoan-mixtures 
#### Mixtures of trees applied to the root of the tree of all animals

Caitlin Cherryh

November 2023

***
### Summary

This github repository contains scripts used to:

1. Estimate trees from 14 empirical phylogenetic datasets with 26 models of substitution
2. Estimate constrained trees for 5 alternate topologies of the Metazoan tree
3. Apply the MAST model and the AU-test to evaluate a multi-tree model for the Metazoan taxa

If you replicate any part of these analyses or use functions from these scripts, please cite this repository.

#### Contents
+ Scripts
     + All scripts necessary to completely replicate this analysis are included in the `code/` folder
    + Each script includes an overview, a list of necessary parameters or file paths,  and a list of software necessary to run that script
+ Output
    + Contains output `.csv` files generated throughout the project
        + `Cherryh_MAST_metazoa_taxa_collation.csv`: table of the name and clade of each species in each dataset
        + `alignment_dimensions.csv`: number of sites and taxa in each alignment    
+ Trees
    + Maximum Likelihood trees
        + Contains all 364 trees generated throughout this process (14 datasets, 26 models of substitution)
    + Hypothesis trees
        + Constrained maximum likelihood trees
+ Taxa reconciliation
    + Table used to make taxa names consistent across datasets
+ Conda enviroment
    + The `environment.yml` file is included to replicate the conda environment used for this project

***
### Instructions to reproduce the analyses:


***
### Treelikeness metrics
