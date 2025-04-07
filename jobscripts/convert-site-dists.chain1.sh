#!/bin/bash
#
#SBATCH --job-name=sitefreqs_chain1
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


cd /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles

python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Borowiec2015/Borowiec2015.Best108.PMSF_Phylobayes_CAT-POISSON.CTEN.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Borowiec2015/Borowiec2015.Best108.PMSF_Phylobayes_CAT-POISSON.PORI.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Chang2015/Chang2015.Chang_AA.PMSF_Phylobayes_CAT-POISSON.CTEN.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Chang2015/Chang2015.Chang_AA.PMSF_Phylobayes_CAT-POISSON.PORI.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Dunn2008/Dunn2008.Dunn2008_FixedNames.PMSF_Phylobayes_CAT-POISSON.CTEN.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Dunn2008/Dunn2008.Dunn2008_FixedNames.PMSF_Phylobayes_CAT-POISSON.PORI.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Laumer2018/Laumer2018.Tplx_phylo_d1.PMSF_Phylobayes_CAT-POISSON.CTEN.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Laumer2018/Laumer2018.Tplx_phylo_d1.PMSF_Phylobayes_CAT-POISSON.PORI.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Laumer2019/Laumer2019.nonbilateria_MARE_BMGE.PMSF_Phylobayes_CAT-POISSON.CTEN.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Laumer2019/Laumer2019.nonbilateria_MARE_BMGE.PMSF_Phylobayes_CAT-POISSON.PORI.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Moroz2014/Moroz2014.ED3d.PMSF_Phylobayes_CAT-POISSON.CTEN.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Moroz2014/Moroz2014.ED3d.PMSF_Phylobayes_CAT-POISSON.PORI.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Nosenko2013/Nosenko2013.nonribosomal_9187_smatrix.PMSF_Phylobayes_CAT-POISSON.CTEN.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Nosenko2013/Nosenko2013.nonribosomal_9187_smatrix.PMSF_Phylobayes_CAT-POISSON.PORI.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Nosenko2013/Nosenko2013.ribosomal_14615_smatrix.PMSF_Phylobayes_CAT-POISSON.CTEN.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Nosenko2013/Nosenko2013.ribosomal_14615_smatrix.PMSF_Phylobayes_CAT-POISSON.PORI.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Philippe2009/Philippe2009.Philippe_etal_superalignment_FixedNames.PMSF_Phylobayes_CAT-POISSON.CTEN.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Philippe2009/Philippe2009.Philippe_etal_superalignment_FixedNames.PMSF_Phylobayes_CAT-POISSON.PORI.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Philippe2011/Philippe2011.UPDUNN_MB_FixedNames.PMSF_Phylobayes_CAT-POISSON.CTEN.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Philippe2011/Philippe2011.UPDUNN_MB_FixedNames.PMSF_Phylobayes_CAT-POISSON.PORI.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Pick2010/Pick2010.Pick2010.PMSF_Phylobayes_CAT-POISSON.CTEN.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Pick2010/Pick2010.Pick2010.PMSF_Phylobayes_CAT-POISSON.PORI.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Ryan2013/Ryan2013.REA_EST_includingXenoturbella.PMSF_Phylobayes_CAT-POISSON.CTEN.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Ryan2013/Ryan2013.REA_EST_includingXenoturbella.PMSF_Phylobayes_CAT-POISSON.PORI.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Whelan2015/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.PMSF_Phylobayes_CAT-POISSON.CTEN.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Whelan2015/Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.PMSF_Phylobayes_CAT-POISSON.PORI.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Whelan2017/Whelan2017.Metazoa_Choano_RCFV_strict.PMSF_Phylobayes_CAT-POISSON.CTEN.chain1.siteprofiles
python3 /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/drenal_cat-pmsf-paper/scripts/convert-site-dists.py /mnt/data/dayhoff/home/u5348329/metazoan-mixtures/CAT_PMSF/01_CAT-PMSF_profiles/Whelan2017/Whelan2017.Metazoa_Choano_RCFV_strict.PMSF_Phylobayes_CAT-POISSON.PORI.chain1.siteprofiles
