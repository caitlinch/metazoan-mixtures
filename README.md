# metazoan-mixtures 
#### Mixtures of trees applied to the root of the tree of all animals

Caitlin Cherryh

November 2023

***
### Summary

This github repository contains scripts used to:

1. Simulate data with different underlying levels of treelikeness and incomplete lineage sorting (ILS)
2. Benchmark existing metrics for treelikeness against these simulations
3. Introduce a new metric for treelikeness in phylogenetic datasets, called the tree proportion

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
