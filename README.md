# metazoan-mixtures 
#### Applying the MAST model to the Animal tree of life

Caitlin Cherryh

November 2023

***
### Summary

This github repository contains scripts used to:

1. Estimate trees from 14 empirical phylogenetic datasets with 26 models of substitution
2. Estimate constrained trees for 5 alternate topologies of the Metazoan tree
3. Apply the [MAST model](https://www.biorxiv.org/content/10.1101/2022.10.06.511210v1) and the AU-test to evaluate a multi-tree model for the Metazoan taxa

If you replicate any part of these analyses or use functions from these scripts, please cite this repository.

***
### Contents
+ Scripts
   + All scripts necessary to completely replicate this analysis are included in the `code/` folder
   + Each script includes an overview, a list of necessary parameters or file paths,  and a list of software necessary to run that script
+ Output
   + Contains output `.csv` files generated throughout the project
      + `Cherryh_MAST_metazoa_taxa_collation.csv`: table of the name and clade of each species in each dataset
      + `alignment_dimensions.csv`: number of sites and taxa in each alignment  
+ Trees
   + Maximum Likelihood trees
      + Contains all 364 trees generated throughout this process (14 datasets, 26 models of substitution)
   + Hypothesis trees
      + Constrained maximum likelihood trees
+ Taxa reconciliation
   + Table used to make taxa names consistent across datasets
+ Conda enviroment
   + The `environment.yml` file is included to replicate the conda environment used for this project

***
### Instructions to reproduce the analyses:
1. Download and install the software programs necessary to repeat these analyses:
   + [IQ-Tree2](http://www.iqtree.org/)
2. Estimate trees
   + Estimate maximum likelihood trees with standard IQ-Tree protein models and profile mixture (PM) models in IQ-Tree using the script `01_estimate_all_trees_parallel.R`
   + Estimate trees with the posterior mean site frequency (PMSF) in IQ-Tree using the script `01_estimate_PMSF_trees.R`
   + To rename tips in all trees to be consistent across datasets, use the script `util_tree_processing.R` 
3. Estimate constrained trees using the best models of evolution in each class using the script `02_estimate_hypothesis_trees.R`
4. Apply the mixture of trees model using the script `03_TreeMixtures.R`
5. Format output csvs using the script `04_reformat_output_dataframes.R`
6. Plot results using the scripts `05_plots.R` and `05_plots_5trees.R`

***
### Datasets
| Original Publication | Repository | Matrix |
| ----------- | ----------- | ----------- |
| Dunn _et al._ (2008) | Li _et al._ (2020) | Dunn2008 |
| Hejnol _et al._ (2009) | Hejnol _et al._ (2009)  | Hejnol_etal_2009 |
| Philippe _et al._ (2009) | Philippe _et al._ (2009) | Philippe_etal_superalignment |
| Pick _et al._ (2010)  | Li _et al._ (2020) | Pick2010 |
| Philippe _et al._ (2011) | Philippe _et al._ (2011) | UPDUNN_MB |
| Nosenko _et al._ (2013a) | Nosenko _et al._ (2013b) | nonribosomal_9187_smatrix |
| Nosenko _et al._ (2013a) | Nosenko _et al._ (2013b) | ribosomal_14615_smatrix |
| Ryan _et al._ (2013)  | Redmond and McLysaght (2021) | REA_alignment_includingXenoturbella |
| Moroz _et al._ (2014) | Li _et al._ (2020) | ED3d |
| Borowiec _et al._ (2015) | Borowiec _et al._ (2016) | Best108 |
| Chang _et al._ (2015) | Feuda _et al._ (2017)  | Chang_AA |
| Whelan _et al._ (2015) | Whelan _et al._ (2016) | Dataset10 |
| Whelan _et al._ (2017ba | Whelan _et al._ (2017b) | Metazoa_Choano_RCFV_strict |
| Laumer _et al._ (2018a) | Laumer _et al._ (2018b) | Tplx_phylo_d1 |
| Laumer _et al._ (2019a) | Laumer _et al._ (2019b) | nonbilateria_MARE_BMGE |

#### Citations
* Borowiec, M.L., Lee, E.K., Chiu, J.C., Plachetzki, D.C., 2016. Data from: Extracting phylogenetic signal and accounting for bias in whole-genome data sets supports the Ctenophora as sister to remaining Metazoa. https://doi.org/10.5061/DRYAD.K6TQ2 
* Borowiec, M.L., Lee, E.K., Chiu, J.C., Plachetzki, D.C., 2015. Extracting phylogenetic signal and accounting for bias in whole-genome data sets supports the Ctenophora as sister to remaining Metazoa. BMC Genomics 16, 987. https://doi.org/10.1186/s12864-015-2146-4 
* Chang, E.S., Neuhof, M., Rubinstein, N.D., Diamant, A., Philippe, H., Huchon, D., Cartwright, P., 2015. Genomic insights into the evolutionary origin of Myxozoa within Cnidaria. Proceedings of the National Academy of Sciences 112, 14912. https://doi.org/10.1073/pnas.1511468112 
* Dunn, C.W., Hejnol, A., Matus, D.Q., Pang, K., Browne, W.E., Smith, S.A., Seaver, E., Rouse, G.W., Obst, M., Edgecombe, G.D., Sørensen, M.V., Haddock, S.H.D., Schmidt-Rhaesa, A., Okusu, A., Kristensen, R.M., Wheeler, W.C., Martindale, M.Q., Giribet, G., 2008. Broad phylogenomic sampling improves resolution of the animal tree of life. Nature 452, 745–749. https://doi.org/10.1038/nature06614 
* Feuda, R., Dohrmann, M., Pett, W., Philippe, H., Rota-Stabelli, O., Lartillot, N., Wörheide, G., Pisani, D., 2017. Data repository for “Improved Modeling of Compositional Heterogeneity Supports Sponges as Sister to All Other Animals.” https://doi.org/10.1016/j.cub.2017.11.008 
* Hejnol, A., Obst, M., Stamatakis, A., Ott, M., Rouse, G.W., Edgecombe, G.D., Martinez, P., Baguñà, J., Bailly, X., Jondelius, U., Wiens, M., Müller, W.E.G., Seaver, E., Wheeler, W.C., Martindale, M.Q., Giribet, G., Dunn, C.W., 2009. Assessing the root of bilaterian animals with scalable phylogenomic methods. Proceedings of the Royal Society B: Biological Sciences 276, 4261–4270. https://doi.org/10.1098/rspb.2009.0896 
* Laumer, C.E., Fernández, R., Lemer, S., Combosch, D., Kocot, K.M., Riesgo, A., Andrade, S.C.S., Sterrer, W., Sørensen, M.V., Giribet, G., 2019a. Revisiting metazoan phylogeny with genomic sampling of all phyla. Proceedings of the Royal Society B: Biological Sciences 286, 20190831. https://doi.org/10.1098/rspb.2019.0831 
* Laumer, C.E., Fernández, R., Lemer, S., Combosch, D., Kocot, K.M., Riesgo, A., Andrade, S.C.S., Sterrer, W., Sørensen, M.V., Giribet, G., 2019b. Data from: Revisiting metazoan phylogeny with genomic sampling of all phyla. https://doi.org/10.5061/DRYAD.293KP3D 
* Laumer, C.E., Gruber-Vodicka, H., Hadfield, M.G., Pearse, V.B., Riesgo, A., Marioni, J.C., Giribet, G., 2018a. Support for a clade of Placozoa and Cnidaria in genes with minimal compositional bias. eLife 7, e36278. https://doi.org/10.7554/eLife.36278 
* Laumer, C.E., Gruber-Vodicka, H., Hadfield, M.G., Pearse, V.B., Riesgo, A., Marioni, J.C., Giribet, G., 2018b. Data from: Support for a clade of Placozoa and Cnidaria in genes with minimal compositional bias. https://doi.org/10.5061/DRYAD.6CM1166 
* Li, Y., Shen, X.-X., Evans, B., Dunn, C.W., Rokas, A., 2020. Data repository for “Rooting the animal tree of life.” https://doi.org/10.6084/m9.figshare.13122881.v1 
* Moroz, L.L., Kocot, K.M., Citarella, M.R., Dosung, S., Norekian, T.P., Povolotskaya, I.S., Grigorenko, A.P., Dailey, C., Berezikov, E., Buckley, K.M., Ptitsyn, A., Reshetov, D., Mukherjee, K., Moroz, T.P., Bobkova, Y., Yu, F., Kapitonov, V.V., Jurka, J., Bobkov, Y.V., Swore, J.J., Girardo, D.O., Fodor, A., Gusev, F., Sanford, R., Bruders, R., Kittler, E., Mills, C.E., Rast, J.P., Derelle, R., Solovyev, V.V., Kondrashov, F.A., Swalla, B.J., Sweedler, J.V., Rogaev, E.I., Halanych, K.M., Kohn, A.B., 2014. The ctenophore genome and the evolutionary origins of neural systems. Nature 510, 109–114. https://doi.org/10.1038/nature13400 
* Nosenko, T., Schreiber, F., Adamska, M., Adamski, M., Eitel, M., Hammel, J., Maldonado, M., Müller, W.E.G., Nickel, M., Schierwater, B., Vacelet, J., Wiens, M., Wörheide, G., 2013a. Deep metazoan phylogeny: when different genes tell different stories. Molecular Phylogenetics and Evolution 67, 223–233. https://doi.org/10.1016/j.ympev.2013.01.010 
* Nosenko, T., Schreiber, F., Adamska, M., Adamski, M., Eitel, M., Hammel, J., Maldonado, M., Müller, W.E.G., Nickel, M., Schierwater, B., Vacelet, J., Wiens, M., Wörheide, G., 2013b. Additional data to: Deep metazoan phylogeny: When different genes tell different stories. https://doi.org/10.5282/ubm/data.55 
* Philippe, H., Brinkmann, H., Lavrov, D.V., Littlewood, D.T.J., Manuel, M., Wörheide, G., Baurain, D., 2011. Resolving difficult phylogenetic questions: why more sequences are not enough. PLoS Biol 9, e1000602–e1000602. https://doi.org/10.1371/journal.pbio.1000602 
* Philippe, H., Derelle, R., Lopez, P., Pick, K., Borchiellini, C., Boury-Esnault, N., Vacelet, J., Renard, E., Houliston, E., Quéinnec, E., Da Silva, C., Wincker, P., Le Guyader, H., Leys, S., Jackson, D.J., Schreiber, F., Erpenbeck, D., Morgenstern, B., Wörheide, G., Manuel, M., 2009. Phylogenomics revives traditional views on deep animal relationships. Current Biology 19, 706–712. https://doi.org/10.1016/j.cub.2009.02.052 
* Pick, K.S., Philippe, H., Schreiber, F., Erpenbeck, D., Jackson, D.J., Wrede, P., Wiens, M., Alié, A., Morgenstern, B., Manuel, M., Wörheide, G., 2010. Improved phylogenomic taxon sampling noticeably affects nonbilaterian relationships. Molecular Biology and Evolution 27, 1983–1987. https://doi.org/10.1093/molbev/msq089 
* Redmond, A.K., McLysaght, A., 2021. FROM: Evidence for sponges as sister to all other animals from partitioned phylogenomics with mixture models and recoding. https://doi.org/10.6084/m9.figshare.12746972.v2 
* Ryan, J.F., Pang, K., Schnitzler, C.E., Nguyen, A.-D., Moreland, R.T., Simmons, D.K., Koch, B.J., Francis, W.R., Havlak, P., Smith, S.A., Putnam, N.H., Haddock, S.H.D., Dunn, C.W., Wolfsberg, T.G., Mullikin, J.C., Martindale, M.Q., Baxevanis, A.D., 2013. The genome of the ctenophore Mnemiopsis leidyi and its implications for cell type evolution. Science 342, 1242592. https://doi.org/10.1126/science.1242592 
* Whelan, N.V., Kocot, K.M., Moroz, L.L., Halanych, K.M., 2016. Error, signal, and the placement of Ctenophora sister to all other animals. v3. https://doi.org/10.6084/m9.figshare.1334306.v3 
* Whelan, N.V., Kocot, K.M., Moroz, L.L., Halanych, K.M., 2015. Error, signal, and the placement of Ctenophora sister to all other animals. Proceedings of the National Academy of Sciences of the United States of America 112, 5773–5778. https://doi.org/10.1073/pnas.1503453112


* Whelan, N.V., Kocot, K.M., Moroz, T.P., Mukherjee, K., Williams, P., Paulay, G., Moroz, L.L., Halanych, K.M., 2017a. Ctenophore relationships and their placement as the sister group to all other animals. Nature Ecology & Evolution 1, 1737–1746. https://doi.org/10.1038/s41559-017-0331-3
* Whelan, N.V., Kocot, K.M., Moroz, T.P., Mukherjee, K., Williams, P., Paulay, G., Moroz, L.L., Halanych, K.M., 2017b. Ctenophora Phylogeny Datasets and Core Orthologs. Dataset. https://doi.org/10.6084/m9.figshare.4484138.v1 
