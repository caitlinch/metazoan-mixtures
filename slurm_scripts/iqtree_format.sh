#!/bin/bash
#
#SBATCH --job-name=tree_mixtures.${id}
#SBATCH --output=/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/%j.${id}.%x.out
#SBATCH --error=/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/%j.${id}.%x.err 
#SBATCH --partition=Standard
#
#SBATCH --time=5:00:00 # 5 hours - needed ~4 for the Whelan2017 alignment
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5000 # 5GB - needed ~4.6GB for the Whelan2017 alignment
#
#SBATCH --mail-user u5348329@anu.edu.au
#SBATCH --mail-type ALL

# Change to working directory for this project
cd /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/


# Run iqtree command
