#!/bin/bash
#
#SBATCH --job-name=tree_mixtures.17
#SBATCH --output=/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/%j.17.%x.out
#SBATCH --error=/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/%j.17.%x.err 
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
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'CAT'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.CAT 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'CF4'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.CF4 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'EHO'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.EHO 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'EX_EHO'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.EX_EHO 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'EX2'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.EX2 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'EX3'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.EX3 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'F81'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.F81 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'GTR'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.GTR 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'GTR20'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.GTR20 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'JTT'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.JTT 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'JTTDCMut'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.JTTDCMut 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'LG'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.LG 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'LG4M'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.LG4M 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'mtZOA'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.mtZOA 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'PMB'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.PMB 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'Poisson'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.Poisson 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'rtREV'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.rtREV 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'UL2'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.UL2 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'UL3'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.UL3 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -mset 'WAG'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.WAG 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa.alignment.phy -m MFP  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.ModelFinder 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2017.Metazoa_Choano_RCFV_strict.aa.alignment.fa -mset 'C10'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2017.Metazoa_Choano_RCFV_strict.C10 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2017.Metazoa_Choano_RCFV_strict.aa.alignment.fa -mset 'C20'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2017.Metazoa_Choano_RCFV_strict.C20 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Whelan2017.Metazoa_Choano_RCFV_strict.aa.alignment.fa -mset 'C30'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Whelan2017.Metazoa_Choano_RCFV_strict.C30 
