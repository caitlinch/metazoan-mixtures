#!/bin/bash
#
#SBATCH --job-name=pb_Pick2010
#SBATCH --output=/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/job_files/%j.%x.out
#SBATCH --error=/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/job_files/%j.%x.err 
#SBATCH --partition=Standard
#
#SBATCH --time=168:00:00 # Nosenko 2013 non-ribo took ~48 hours to run ~35 generations
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=5000 # 5GB - Nosenko 2013 non-ribo used 1.12GB for 1 chain of PhyloBayes
#
#SBATCH --mail-user u5348329@anu.edu.au
#SBATCH --mail-type TIME_LIMIT,FAIL

# Make dataset directory
mkdir -p /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Pick2010/

# Run PhyloBayes pb
/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/phylobayes/phylobayes-4.1e/data/pb -s /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -T /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/hypothesis_trees/CXX/Pick2010.Pick2010.LG_C60.ML_H1.treefile -cat -poisson Pick2010.Pick2010.PMSF_Phylobayes_CAT-POISSON.CTEN.chain1
/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/phylobayes/phylobayes-4.1e/data/pb -s /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -T /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/hypothesis_trees/CXX/Pick2010.Pick2010.LG_C60.ML_H1.treefile -cat -poisson Pick2010.Pick2010.PMSF_Phylobayes_CAT-POISSON.CTEN.chain2

# Run PhyloBayes tracecomp
/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/phylobayes/phylobayes-4.1e/data/tracecomp -x 500 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Pick2010/Pick2010.Pick2010.PMSF_Phylobayes_CAT-POISSON.CTEN.chain1 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Pick2010/Pick2010.Pick2010.PMSF_Phylobayes_CAT-POISSON.CTEN.chain2

# Run PhyloBayes readpb
/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/phylobayes/phylobayes-4.1e/data/readpb -ss -x 500 10 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Pick2010/Pick2010.Pick2010.PMSF_Phylobayes_CAT-POISSON.CTEN.chain1
/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/phylobayes/phylobayes-4.1e/data/readpb -ss -x 500 10 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Pick2010/Pick2010.Pick2010.PMSF_Phylobayes_CAT-POISSON.CTEN.chain2

