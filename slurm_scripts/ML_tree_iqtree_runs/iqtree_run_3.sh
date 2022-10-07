#!/bin/bash
#
#SBATCH --job-name=tree_mixtures.3
#SBATCH --output=/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/%j.3.%x.out
#SBATCH --error=/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/%j.3.%x.err 
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
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Chang2015.Chang_AA.aa.alignment.phy -mset 'Poisson'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Chang2015.Chang_AA.Poisson 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Chang2015.Chang_AA.aa.alignment.phy -mset 'rtREV'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Chang2015.Chang_AA.rtREV 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Chang2015.Chang_AA.aa.alignment.phy -mset 'UL2'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Chang2015.Chang_AA.UL2 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Chang2015.Chang_AA.aa.alignment.phy -mset 'UL3'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Chang2015.Chang_AA.UL3 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Chang2015.Chang_AA.aa.alignment.phy -mset 'WAG'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Chang2015.Chang_AA.WAG 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Chang2015.Chang_AA.aa.alignment.phy -m MFP  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Chang2015.Chang_AA.ModelFinder 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Dunn2008.Dunn2008_FixedNames.aa.alignment.fasta -mset 'C10'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Dunn2008.Dunn2008_FixedNames.C10 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Dunn2008.Dunn2008_FixedNames.aa.alignment.fasta -mset 'C20'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Dunn2008.Dunn2008_FixedNames.C20 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Dunn2008.Dunn2008_FixedNames.aa.alignment.fasta -mset 'C30'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Dunn2008.Dunn2008_FixedNames.C30 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Dunn2008.Dunn2008_FixedNames.aa.alignment.fasta -mset 'C40'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Dunn2008.Dunn2008_FixedNames.C40 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Dunn2008.Dunn2008_FixedNames.aa.alignment.fasta -mset 'C50'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Dunn2008.Dunn2008_FixedNames.C50 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Dunn2008.Dunn2008_FixedNames.aa.alignment.fasta -mset 'C60'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Dunn2008.Dunn2008_FixedNames.C60 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Dunn2008.Dunn2008_FixedNames.aa.alignment.fasta -mset 'CAT'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Dunn2008.Dunn2008_FixedNames.CAT 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Dunn2008.Dunn2008_FixedNames.aa.alignment.fasta -mset 'CF4'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Dunn2008.Dunn2008_FixedNames.CF4 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Dunn2008.Dunn2008_FixedNames.aa.alignment.fasta -mset 'EHO'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Dunn2008.Dunn2008_FixedNames.EHO 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Dunn2008.Dunn2008_FixedNames.aa.alignment.fasta -mset 'EX_EHO'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Dunn2008.Dunn2008_FixedNames.EX_EHO 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Dunn2008.Dunn2008_FixedNames.aa.alignment.fasta -mset 'EX2'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Dunn2008.Dunn2008_FixedNames.EX2 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Dunn2008.Dunn2008_FixedNames.aa.alignment.fasta -mset 'EX3'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Dunn2008.Dunn2008_FixedNames.EX3 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Dunn2008.Dunn2008_FixedNames.aa.alignment.fasta -mset 'F81'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Dunn2008.Dunn2008_FixedNames.F81 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Dunn2008.Dunn2008_FixedNames.aa.alignment.fasta -mset 'GTR'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Dunn2008.Dunn2008_FixedNames.GTR 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Dunn2008.Dunn2008_FixedNames.aa.alignment.fasta -mset 'GTR20'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Dunn2008.Dunn2008_FixedNames.GTR20 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Dunn2008.Dunn2008_FixedNames.aa.alignment.fasta -mset 'JTT'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Dunn2008.Dunn2008_FixedNames.JTT 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Dunn2008.Dunn2008_FixedNames.aa.alignment.fasta -mset 'JTTDCMut'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Dunn2008.Dunn2008_FixedNames.JTTDCMut 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Dunn2008.Dunn2008_FixedNames.aa.alignment.fasta -mset 'LG'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Dunn2008.Dunn2008_FixedNames.LG 
