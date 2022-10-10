#!/bin/bash
#
#SBATCH --job-name=tree_mixtures.${id}
#SBATCH --output=/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/%j.${id}.%x.out
#SBATCH --error=/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/%j.${id}.%x.err 
#SBATCH --partition=Standard
#
#SBATCH --time=72:00:00 # 72 hours = 3 days
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=7000 # 7GB (largest alignment run so far needed 6.5 GB)
#
#SBATCH --mail-user u5348329@anu.edu.au
#SBATCH --mail-type NONE

# Change to working directory for this project
cd /home/u5348329/metazoan-mixtures/output/maximum_likelihood_trees/

# Run iqtree command
