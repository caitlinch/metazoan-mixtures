#!/bin/bash
#
#SBATCH --job-name=tree_mixtures.14
#SBATCH --output=/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/%j.14.%x.out
#SBATCH --error=/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/%j.14.%x.err 
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
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'JTT'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.JTT 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'JTTDCMut'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.JTTDCMut 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'LG'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.LG 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'LG4M'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.LG4M 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'mtZOA'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.mtZOA 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'PMB'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.PMB 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'Poisson'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.Poisson 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'rtREV'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.rtREV 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'UL2'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.UL2 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'UL3'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.UL3 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'WAG'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.WAG 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -m MFP  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.ModelFinder 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Ryan2013.REA_EST_includingXenoturbella.aa.alignment.fas -mset 'C10'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Ryan2013.REA_EST_includingXenoturbella.C10 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Ryan2013.REA_EST_includingXenoturbella.aa.alignment.fas -mset 'C20'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Ryan2013.REA_EST_includingXenoturbella.C20 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Ryan2013.REA_EST_includingXenoturbella.aa.alignment.fas -mset 'C30'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Ryan2013.REA_EST_includingXenoturbella.C30 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Ryan2013.REA_EST_includingXenoturbella.aa.alignment.fas -mset 'C40'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Ryan2013.REA_EST_includingXenoturbella.C40 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Ryan2013.REA_EST_includingXenoturbella.aa.alignment.fas -mset 'C50'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Ryan2013.REA_EST_includingXenoturbella.C50 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Ryan2013.REA_EST_includingXenoturbella.aa.alignment.fas -mset 'C60'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Ryan2013.REA_EST_includingXenoturbella.C60
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Ryan2013.REA_EST_includingXenoturbella.aa.alignment.fas -mset 'CF4'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Ryan2013.REA_EST_includingXenoturbella.CF4 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Ryan2013.REA_EST_includingXenoturbella.aa.alignment.fas -mset 'EHO'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Ryan2013.REA_EST_includingXenoturbella.EHO 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Ryan2013.REA_EST_includingXenoturbella.aa.alignment.fas -mset 'EX_EHO'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Ryan2013.REA_EST_includingXenoturbella.EX_EHO 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Ryan2013.REA_EST_includingXenoturbella.aa.alignment.fas -mset 'EX2'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Ryan2013.REA_EST_includingXenoturbella.EX2 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Ryan2013.REA_EST_includingXenoturbella.aa.alignment.fas -mset 'EX3'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Ryan2013.REA_EST_includingXenoturbella.EX3 
