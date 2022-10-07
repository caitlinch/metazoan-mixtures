#!/bin/bash
#
#SBATCH --job-name=tree_mixtures
#SBATCH --output=%j.%x.out
#SBATCH --error=%j.%x.err 
#SBATCH --partition=Standard
#
#SBATCH --time=3:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=100
#
##SBATCH --mail-user u5348329@anu.edu.au
##SBATCH --mail-type ALL

# Change to working directory for this project
cd /home/u5348329/metazoan-mixtures/

# Activate Anaconda work environment
source /home/u5348329/.bashrc
source conda activate tree_mixtures

# Run Rscript in a clean R instance
Rscript --vanilla --verbose --no-restore --no-save code/01_estimate_all_trees_parallel.R &> tree_mixtures_R.log
