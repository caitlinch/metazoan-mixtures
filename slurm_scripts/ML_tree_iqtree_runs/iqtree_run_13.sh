#!/bin/bash
#
#SBATCH --job-name=tree_mixtures.13
#SBATCH --output=/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/%j.13.%x.out
#SBATCH --error=/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/%j.13.%x.err 
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
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Philippe2011.UPDUNN_MB_FixedNames.aa.alignment.nex -mset 'LG4M'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Philippe2011.UPDUNN_MB_FixedNames.LG4M 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Philippe2011.UPDUNN_MB_FixedNames.aa.alignment.nex -mset 'mtZOA'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Philippe2011.UPDUNN_MB_FixedNames.mtZOA 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Philippe2011.UPDUNN_MB_FixedNames.aa.alignment.nex -mset 'PMB'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Philippe2011.UPDUNN_MB_FixedNames.PMB 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Philippe2011.UPDUNN_MB_FixedNames.aa.alignment.nex -mset 'Poisson'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Philippe2011.UPDUNN_MB_FixedNames.Poisson 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Philippe2011.UPDUNN_MB_FixedNames.aa.alignment.nex -mset 'rtREV'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Philippe2011.UPDUNN_MB_FixedNames.rtREV 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Philippe2011.UPDUNN_MB_FixedNames.aa.alignment.nex -mset 'UL2'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Philippe2011.UPDUNN_MB_FixedNames.UL2 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Philippe2011.UPDUNN_MB_FixedNames.aa.alignment.nex -mset 'UL3'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Philippe2011.UPDUNN_MB_FixedNames.UL3 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Philippe2011.UPDUNN_MB_FixedNames.aa.alignment.nex -mset 'WAG'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Philippe2011.UPDUNN_MB_FixedNames.WAG 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Philippe2011.UPDUNN_MB_FixedNames.aa.alignment.nex -m MFP  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Philippe2011.UPDUNN_MB_FixedNames.ModelFinder 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'C10'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.C10 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'C20'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.C20 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'C30'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.C30 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'C40'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.C40 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'C50'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.C50 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'C60'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.C60 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'CAT'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.CAT 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'CF4'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.CF4 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'EHO'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.EHO 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'EX_EHO'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.EX_EHO 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'EX2'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.EX2 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'EX3'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.EX3 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'F81'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.F81 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'GTR'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.GTR 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Pick2010.Pick2010.aa.alignment.phy -mset 'GTR20'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Pick2010.Pick2010.GTR20 
