#!/bin/bash
#
#SBATCH --job-name=Ryan2013_c1
#SBATCH --output=/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/job_files/%j.%x.out
#SBATCH --error=/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/job_files/%j.%x.err 
#SBATCH --partition=Standard
#
#SBATCH --time=336:00:00 # Nosenko 2013 non-ribo took ~48 hours to run ~35 generations - give this dataset 2 weeks
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=5000 # 5GB - Nosenko 2013 non-ribo used 1.12GB for 1 chain of PhyloBayes
#
#SBATCH --mail-user u5348329@anu.edu.au
#SBATCH --mail-type TIME_LIMIT,FAIL

# Change to directory
cd /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/

# Run PhyloBayes pb
# /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/phylobayes/phylobayes-4.1e/data/pb \
#    -d /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/data_all/Ryan2013.REA_EST_includingXenoturbella.aa.alignment.fas \
#    -T /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/hypothesis_trees/CXX/Ryan2013.REA_EST_includingXenoturbella.LG_C60.ML_H1.treefile \
#    -cat -poisson \
#    Ryan2013.REA_EST_includingXenoturbella.PMSF_Phylobayes_CAT-POISSON.CTEN.100K.chain1
/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/phylobayes/phylobayes-4.1e/data/pb \
    -d /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/data_all/Ryan2013.REA_EST_includingXenoturbella.aa.alignment.fas \
    -T /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/hypothesis_trees/CXX/Ryan2013.REA_EST_includingXenoturbella.LG_C60.ML_H1.treefile \
    -cat -poisson \
    Ryan2013.REA_EST_includingXenoturbella.PMSF_Phylobayes_CAT-POISSON.CTEN.100K.chain2

# Run PhyloBayes tracecomp
# /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/phylobayes/phylobayes-4.1e/data/tracecomp \
#    -x 500 \
#    /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Ryan2013/Ryan2013.REA_EST_includingXenoturbella.PMSF_Phylobayes_CAT-POISSON.CTEN.100K.chain1 \
#    /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Ryan2013/Ryan2013.REA_EST_includingXenoturbella.PMSF_Phylobayes_CAT-POISSON.CTEN.100K.chain2

# Run PhyloBayes readpb
# /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/phylobayes/phylobayes-4.1e/data/readpb \
#    -ss -x 500 10 \
#    /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Ryan2013/Ryan2013.REA_EST_includingXenoturbella.PMSF_Phylobayes_CAT-POISSON.CTEN.100K.chain1
# /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/phylobayes/phylobayes-4.1e/data/readpb \
#    -ss -x 500 10 \
#    /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Ryan2013/Ryan2013.REA_EST_includingXenoturbella.PMSF_Phylobayes_CAT-POISSON.CTEN.100K.chain2

