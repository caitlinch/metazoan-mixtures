#!/bin/bash
#
#SBATCH --job-name=tree_mixtures.6
#SBATCH --output=/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/%j.6.%x.out
#SBATCH --error=/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/%j.6.%x.err 
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
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Laumer2018.Tplx_phylo_d1.aa.alignment.phylip -mset 'F81'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Laumer2018.Tplx_phylo_d1.F81 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Laumer2018.Tplx_phylo_d1.aa.alignment.phylip -mset 'GTR'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Laumer2018.Tplx_phylo_d1.GTR 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Laumer2018.Tplx_phylo_d1.aa.alignment.phylip -mset 'GTR20'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Laumer2018.Tplx_phylo_d1.GTR20 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Laumer2018.Tplx_phylo_d1.aa.alignment.phylip -mset 'JTT'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Laumer2018.Tplx_phylo_d1.JTT 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Laumer2018.Tplx_phylo_d1.aa.alignment.phylip -mset 'JTTDCMut'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Laumer2018.Tplx_phylo_d1.JTTDCMut 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Laumer2018.Tplx_phylo_d1.aa.alignment.phylip -mset 'LG'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Laumer2018.Tplx_phylo_d1.LG 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Laumer2018.Tplx_phylo_d1.aa.alignment.phylip -mset 'LG4M'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Laumer2018.Tplx_phylo_d1.LG4M 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Laumer2018.Tplx_phylo_d1.aa.alignment.phylip -mset 'mtZOA'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Laumer2018.Tplx_phylo_d1.mtZOA 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Laumer2018.Tplx_phylo_d1.aa.alignment.phylip -mset 'PMB'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Laumer2018.Tplx_phylo_d1.PMB 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Laumer2018.Tplx_phylo_d1.aa.alignment.phylip -mset 'Poisson'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Laumer2018.Tplx_phylo_d1.Poisson 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Laumer2018.Tplx_phylo_d1.aa.alignment.phylip -mset 'rtREV'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Laumer2018.Tplx_phylo_d1.rtREV 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Laumer2018.Tplx_phylo_d1.aa.alignment.phylip -mset 'UL2'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Laumer2018.Tplx_phylo_d1.UL2 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Laumer2018.Tplx_phylo_d1.aa.alignment.phylip -mset 'UL3'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Laumer2018.Tplx_phylo_d1.UL3 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Laumer2018.Tplx_phylo_d1.aa.alignment.phylip -mset 'WAG'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Laumer2018.Tplx_phylo_d1.WAG 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Laumer2018.Tplx_phylo_d1.aa.alignment.phylip -m MFP  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Laumer2018.Tplx_phylo_d1.ModelFinder 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Laumer2019.nonbilateria_MARE_BMGE.aa.alignment.phylip -mset 'C10'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Laumer2019.nonbilateria_MARE_BMGE.C10 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Laumer2019.nonbilateria_MARE_BMGE.aa.alignment.phylip -mset 'C20'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Laumer2019.nonbilateria_MARE_BMGE.C20 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Laumer2019.nonbilateria_MARE_BMGE.aa.alignment.phylip -mset 'C30'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Laumer2019.nonbilateria_MARE_BMGE.C30 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Laumer2019.nonbilateria_MARE_BMGE.aa.alignment.phylip -mset 'C40'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Laumer2019.nonbilateria_MARE_BMGE.C40 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Laumer2019.nonbilateria_MARE_BMGE.aa.alignment.phylip -mset 'C50'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Laumer2019.nonbilateria_MARE_BMGE.C50 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Laumer2019.nonbilateria_MARE_BMGE.aa.alignment.phylip -mset 'C60'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Laumer2019.nonbilateria_MARE_BMGE.C60 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Laumer2019.nonbilateria_MARE_BMGE.aa.alignment.phylip -mset 'CAT'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Laumer2019.nonbilateria_MARE_BMGE.CAT 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Laumer2019.nonbilateria_MARE_BMGE.aa.alignment.phylip -mset 'CF4'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Laumer2019.nonbilateria_MARE_BMGE.CF4 
/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2 -s /home/u5348329/metazoan-mixtures/data_all/Laumer2019.nonbilateria_MARE_BMGE.aa.alignment.phylip -mset 'EHO'  -mrate 'E,I,G,I+G,R,I+R'  -bb 1000  -nt 1 -pre Laumer2019.nonbilateria_MARE_BMGE.EHO 
