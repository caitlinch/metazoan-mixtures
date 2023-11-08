# metazoan-mixtures 
#### Mixtures of trees applied to the root of the tree of all animals

Caitlin Cherryh

November 2023

***
### Summary

This github repository contains scripts used to:

1. Estimate trees from 14 empirical phylogenetic datasets with 26 models of substitution
2. Estimate constrained trees for 5 alternate topologies of the Metazoan tree
3. Apply the [MAST model](https://www.biorxiv.org/content/10.1101/2022.10.06.511210v1) and the AU-test to evaluate a multi-tree model for the Metazoan taxa

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
1. Download and install the software programs necessary to repeat these analyses:
    + [IQ-Tree2](http://www.iqtree.org/)
2. Estimate trees
    a. Estimate maximum likelihood trees with standard IQ-Tree protein models and profile mixture (PM) models in IQ-Tree using the script `01_estimate_all_trees_parallel.R`
   b. Estimate trees with the posterior mean site frequency (PMSF) in IQ-Tree using the script `01_estimate_PMSF_trees.R`
   c. To rename tips in all trees to be consistent across datasets, use the script `util_tree_processing.R`   
4. Estimate constrained trees using the best models of evolution in each class using the script `02_estimate_hypothesis_trees.R`
5. Apply the mixture of trees model using the script `03_TreeMixtures.R`
6. Format output csvs using the script `04_reformat_output_dataframes.R`
7. Plot results using the scripts `05_plots.R` and `05_plots_5trees.R`

***
### Datasets
