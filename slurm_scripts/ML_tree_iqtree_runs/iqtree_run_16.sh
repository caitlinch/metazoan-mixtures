#!/bin/bash
#
#SBATCH --job-name=tree_mixtures.16
#SBATCH --output=/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/%j.16.%x.out
#SBATCH --error=/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/%j.16.%x.err 
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
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Simion2017.supermatrix_97sp_401632pos_1719genes.aa.alignment.fasta -mset 'EX_EHO'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Simion2017.supermatrix_97sp_401632pos_1719genes.EX_EHO 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Simion2017.supermatrix_97sp_401632pos_1719genes.aa.alignment.fasta -mset 'EX2'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Simion2017.supermatrix_97sp_401632pos_1719genes.EX2 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Simion2017.supermatrix_97sp_401632pos_1719genes.aa.alignment.fasta -mset 'EX3'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Simion2017.supermatrix_97sp_401632pos_1719genes.EX3 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Simion2017.supermatrix_97sp_401632pos_1719genes.aa.alignment.fasta -mset 'F81'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Simion2017.supermatrix_97sp_401632pos_1719genes.F81 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Simion2017.supermatrix_97sp_401632pos_1719genes.aa.alignment.fasta -mset 'GTR'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Simion2017.supermatrix_97sp_401632pos_1719genes.GTR 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Simion2017.supermatrix_97sp_401632pos_1719genes.aa.alignment.fasta -mset 'GTR20'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Simion2017.supermatrix_97sp_401632pos_1719genes.GTR20 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Simion2017.supermatrix_97sp_401632pos_1719genes.aa.alignment.fasta -mset 'JTT'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Simion2017.supermatrix_97sp_401632pos_1719genes.JTT 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Simion2017.supermatrix_97sp_401632pos_1719genes.aa.alignment.fasta -mset 'JTTDCMut'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Simion2017.supermatrix_97sp_401632pos_1719genes.JTTDCMut 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Simion2017.supermatrix_97sp_401632pos_1719genes.aa.alignment.fasta -mset 'LG'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Simion2017.supermatrix_97sp_401632pos_1719genes.LG 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Simion2017.supermatrix_97sp_401632pos_1719genes.aa.alignment.fasta -mset 'LG4M'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Simion2017.supermatrix_97sp_401632pos_1719genes.LG4M 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Simion2017.supermatrix_97sp_401632pos_1719genes.aa.alignment.fasta -mset 'mtZOA'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Simion2017.supermatrix_97sp_401632pos_1719genes.mtZOA 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Simion2017.supermatrix_97sp_401632pos_1719genes.aa.alignment.fasta -mset 'PMB'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Simion2017.supermatrix_97sp_401632pos_1719genes.PMB 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Simion2017.supermatrix_97sp_401632pos_1719genes.aa.alignment.fasta -mset 'Poisson'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Simion2017.supermatrix_97sp_401632pos_1719genes.Poisson 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Simion2017.supermatrix_97sp_401632pos_1719genes.aa.alignment.fasta -mset 'rtREV'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Simion2017.supermatrix_97sp_401632pos_1719genes.rtREV 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Simion2017.supermatrix_97sp_401632pos_1719genes.aa.alignment.fasta -mset 'UL2'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Simion2017.supermatrix_97sp_401632pos_1719genes.UL2 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Simion2017.supermatrix_97sp_401632pos_1719genes.aa.alignment.fasta -mset 'UL3'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Simion2017.supermatrix_97sp_401632pos_1719genes.UL3 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Simion2017.supermatrix_97sp_401632pos_1719genes.aa.alignment.fasta -mset 'WAG'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Simion2017.supermatrix_97sp_401632pos_1719genes.WAG 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Simion2017.supermatrix_97sp_401632pos_1719genes.aa.alignment.fasta -m MFP  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Simion2017.supermatrix_97sp_401632pos_1719genes.ModelFinder 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'C10'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.C10 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'C20'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.C20 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'C30'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.C30 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'C40'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.C40 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'C50'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.C50 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'C60'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.C60 
