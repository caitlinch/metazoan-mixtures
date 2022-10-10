#!/bin/bash
#
#SBATCH --job-name=tree_mixtures.9
#SBATCH --output=/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/%j.9.%x.out
#SBATCH --error=/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/%j.9.%x.err 
#SBATCH --partition=Standard
#
#SBATCH --time=72:00:00 # 5 hours - needed ~4 for the Whelan2017 alignment (give 72 hours)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5000 # 5GB - needed ~4.6GB for the Whelan2017 alignment
#
#SBATCH --mail-user u5348329@anu.edu.au
#SBATCH --mail-type NONE

# Change to working directory for this project
cd /home/u5348329/metazoan-mixtures/output/maximum_likelihood_trees/

# Run iqtree command
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Nosenko2013.nonribosomal_9187_smatrix.aa.alignment.phy -mset 'C40'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Nosenko2013.nonribosomal_9187_smatrix.C40 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Nosenko2013.nonribosomal_9187_smatrix.aa.alignment.phy -mset 'C50'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Nosenko2013.nonribosomal_9187_smatrix.C50 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Nosenko2013.nonribosomal_9187_smatrix.aa.alignment.phy -mset 'C60'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Nosenko2013.nonribosomal_9187_smatrix.C60
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Nosenko2013.nonribosomal_9187_smatrix.aa.alignment.phy -mset 'CF4'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Nosenko2013.nonribosomal_9187_smatrix.CF4 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Nosenko2013.nonribosomal_9187_smatrix.aa.alignment.phy -mset 'EHO'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Nosenko2013.nonribosomal_9187_smatrix.EHO 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Nosenko2013.nonribosomal_9187_smatrix.aa.alignment.phy -mset 'EX_EHO'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Nosenko2013.nonribosomal_9187_smatrix.EX_EHO 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Nosenko2013.nonribosomal_9187_smatrix.aa.alignment.phy -mset 'EX2'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Nosenko2013.nonribosomal_9187_smatrix.EX2 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Nosenko2013.nonribosomal_9187_smatrix.aa.alignment.phy -mset 'EX3'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Nosenko2013.nonribosomal_9187_smatrix.EX3
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Nosenko2013.nonribosomal_9187_smatrix.aa.alignment.phy -mset 'GTR20'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Nosenko2013.nonribosomal_9187_smatrix.GTR20 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Nosenko2013.nonribosomal_9187_smatrix.aa.alignment.phy -mset 'JTT'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Nosenko2013.nonribosomal_9187_smatrix.JTT 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Nosenko2013.nonribosomal_9187_smatrix.aa.alignment.phy -mset 'JTTDCMut'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Nosenko2013.nonribosomal_9187_smatrix.JTTDCMut 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Nosenko2013.nonribosomal_9187_smatrix.aa.alignment.phy -mset 'LG'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Nosenko2013.nonribosomal_9187_smatrix.LG 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Nosenko2013.nonribosomal_9187_smatrix.aa.alignment.phy -mset 'LG4M'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Nosenko2013.nonribosomal_9187_smatrix.LG4M 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Nosenko2013.nonribosomal_9187_smatrix.aa.alignment.phy -mset 'mtZOA'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Nosenko2013.nonribosomal_9187_smatrix.mtZOA 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Nosenko2013.nonribosomal_9187_smatrix.aa.alignment.phy -mset 'PMB'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Nosenko2013.nonribosomal_9187_smatrix.PMB 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Nosenko2013.nonribosomal_9187_smatrix.aa.alignment.phy -mset 'Poisson'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Nosenko2013.nonribosomal_9187_smatrix.Poisson 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Nosenko2013.nonribosomal_9187_smatrix.aa.alignment.phy -mset 'rtREV'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Nosenko2013.nonribosomal_9187_smatrix.rtREV 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Nosenko2013.nonribosomal_9187_smatrix.aa.alignment.phy -mset 'UL2'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Nosenko2013.nonribosomal_9187_smatrix.UL2 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Nosenko2013.nonribosomal_9187_smatrix.aa.alignment.phy -mset 'UL3'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Nosenko2013.nonribosomal_9187_smatrix.UL3 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Nosenko2013.nonribosomal_9187_smatrix.aa.alignment.phy -mset 'WAG'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Nosenko2013.nonribosomal_9187_smatrix.WAG 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Nosenko2013.nonribosomal_9187_smatrix.aa.alignment.phy -m MFP  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Nosenko2013.nonribosomal_9187_smatrix.ModelFinder 
