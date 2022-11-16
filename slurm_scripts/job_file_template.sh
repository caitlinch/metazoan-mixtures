#!/bin/bash
#
#SBATCH --job-name=tm.${id}
#SBATCH --output=/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/job_files/%j.%x.out
#SBATCH --error=/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/job_files/%j.%x.err 
#SBATCH --partition=Standard
#
#SBATCH --time=120:00:00 # 5 days
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=7000 # 7GB
#
#SBATCH --mail-user u5348329@anu.edu.au
#SBATCH --mail-type TIME_LIMIT,FAIL

# Change to working directory for this project
cd /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/output/maximum_likelihood_trees/

# Run iqtree command
