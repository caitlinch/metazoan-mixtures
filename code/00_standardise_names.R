# metazoan-mixtures/code/00_standardise_names.R
## This script standardises the names for all species for all datasets
# Caitlin Cherryh, 2022

# For this script, you will need the naming csv from Li et. al. (2021)
#     Available from the data repository: https://figshare.com/articles/dataset/Rooting_the_animal_tree_of_life/13122881
#     Download the reconciliation.tar.xz and identify the location of the "reconciliation/taxonomy_info/taxon_table.tsv" file



#### 1. Input parameters ####
## Specify parameters:
# output_dir            <- Directory for IQ-Tree output (trees and tree mixtures)
# repo_dir              <- Location of caitlinch/metazoan-mixtures github repository
# taxon_table_path      <- Location of the reconciliation/taxonomy_info/taxon_table.tsv from Li et. al. (2021)
# manual_taxonomy_path  <- Location of the reconciliation/taxonomy_info/manual_taxonomy_map.tsv from Li et. al. (2021)

location = "local"
if (location == "local"){
  output_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/"
  repo_dir <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"
  taxon_table_path <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/00_Li2021_supp/reconciliation_keep/taxonomy_info/taxon_table.tsv"
  manual_taxonomy_path <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/00_Li2021_supp/reconciliation_keep/taxonomy_info/manual_taxonomy_map.tsv"
} 



#### 2. Open packages and source functions ####
# Source functions and taxa lists
source(paste0(repo_dir, "code/func_naming.R"))
source(paste0(repo_dir, "code/data_dataset_info.R"))

# Remove the individual dataset lists (only need collated lists) (yes this is a bit cheeky)
rm(borowiec2015_list, chang2015_list, dunn2008_list, hejnol2009_list, laumer2018_list, laumer2019_list, moroz2014_list, nosenko2013_list, philippe2009_list,
     philippe2011_list, pick2010_list, ryan2013_list, simion2017_list, whelan2015_list, whelan2017_list, all_models)

#### 3. Prepare name csv ####
# Set a file path for the name csv
mastmet_file_path <- paste0(output_dir, "MAST_metazoa_taxa_collation.csv")
if (file.exists(mastmet_file_path) == FALSE){
  # Create a new data frame with all the taxa from all the matrices you're using
  mastmet_df <- data.frame("dataset" = c(rep("Borowiec0215", extract.taxa.vector(all_datasets[["Borowiec2015"]])$number),
                                         rep("Chang2015", extract.taxa.vector(all_datasets[["Chang2015"]])$number),
                                         rep("Dunn2008", extract.taxa.vector(all_datasets[["Dunn2008"]])$number),
                                         rep("Hejnol2009", extract.taxa.vector(all_datasets[["Hejnol2009"]])$number),
                                         rep("Laumer2018", extract.taxa.vector(filter.matrix.names(all_datasets[["Laumer2018"]], matrix_taxa[["Laumer2018.Tplx_phylo_d1.aa"]]))$number),
                                         rep("Laumer2019", extract.taxa.vector(filter.matrix.names(all_datasets[["Laumer2019"]], matrix_taxa[["Laumer2019.nonbilateria_MARE_BMGE.aa"]]))$number),
                                         rep("Moroz2014", extract.taxa.vector(filter.matrix.names(all_datasets[["Moroz2014"]], matrix_taxa[["Moroz2014.ED3d.aa"]]))$number),
                                         rep("Nosenko2013.nonribosomal", extract.taxa.vector(filter.matrix.names(all_datasets[["Nosenko2013"]], matrix_taxa[["Nosenko2013.nonribosomal_9187_smatrix.aa"]]))$number),
                                         rep("Nosenko2013.ribosomal", extract.taxa.vector(filter.matrix.names(all_datasets[["Nosenko2013"]], matrix_taxa[["Nosenko2013.ribosomal_14615_smatrix.aa"]]))$number),
                                         rep("Philippe2009", extract.taxa.vector(all_datasets[["Philippe2009"]])$number),
                                         rep("Philippe2011", extract.taxa.vector(all_datasets[["Philippe2011"]])$number),
                                         rep("Pick2010", extract.taxa.vector(all_datasets[["Pick2010"]])$number),
                                         rep("Ryan2013", extract.taxa.vector(all_datasets[["Ryan2013"]])$number),
                                         rep("Simion2017", extract.taxa.vector(filter.matrix.names(all_datasets[["Simion2017"]], matrix_taxa[["Simion2017.supermatrix_97sp_401632pos_1719genes.aa"]]))$number),
                                         rep("Whelan2015", extract.taxa.vector(filter.matrix.names(all_datasets[["Whelan2015"]], matrix_taxa[["Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa"]]))$number),
                                         rep("Whelan2017", extract.taxa.vector(all_datasets[["Whelan2017"]])$number)), 
                           "original_name" = c(extract.taxa.vector(all_datasets[["Borowiec2015"]])$taxa, 
                                               extract.taxa.vector(all_datasets[["Chang2015"]])$taxa,
                                               extract.taxa.vector(all_datasets[["Dunn2008"]])$taxa,
                                               extract.taxa.vector(all_datasets[["Hejnol2009"]])$taxa,
                                               extract.taxa.vector(filter.matrix.names(all_datasets[["Laumer2018"]], matrix_taxa[["Laumer2018.Tplx_phylo_d1.aa"]]))$taxa,
                                               extract.taxa.vector(filter.matrix.names(all_datasets[["Laumer2019"]], matrix_taxa[["Laumer2019.nonbilateria_MARE_BMGE.aa"]]))$taxa,
                                               extract.taxa.vector(filter.matrix.names(all_datasets[["Moroz2014"]], matrix_taxa[["Moroz2014.ED3d.aa"]]))$taxa,
                                               extract.taxa.vector(filter.matrix.names(all_datasets[["Nosenko2013"]], matrix_taxa[["Nosenko2013.nonribosomal_9187_smatrix.aa"]]))$taxa,
                                               extract.taxa.vector(filter.matrix.names(all_datasets[["Nosenko2013"]], matrix_taxa[["Nosenko2013.ribosomal_14615_smatrix.aa"]]))$taxa,
                                               extract.taxa.vector(all_datasets[["Philippe2009"]])$taxa,
                                               extract.taxa.vector(all_datasets[["Philippe2011"]])$taxa,
                                               extract.taxa.vector(all_datasets[["Pick2010"]])$taxa,
                                               extract.taxa.vector(all_datasets[["Ryan2013"]])$taxa,
                                               extract.taxa.vector(filter.matrix.names(all_datasets[["Simion2017"]], matrix_taxa[["Simion2017.supermatrix_97sp_401632pos_1719genes.aa"]]))$taxa,
                                               extract.taxa.vector(filter.matrix.names(all_datasets[["Whelan2015"]], matrix_taxa[["Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa"]]))$taxa,
                                               extract.taxa.vector(all_datasets[["Whelan2017"]])$taxa),
                           "clade" = c(extract.taxa.vector(all_datasets[["Borowiec2015"]])$clade, 
                                       extract.taxa.vector(all_datasets[["Chang2015"]])$clade,
                                       extract.taxa.vector(all_datasets[["Dunn2008"]])$clade,
                                       extract.taxa.vector(all_datasets[["Hejnol2009"]])$clade,
                                       extract.taxa.vector(filter.matrix.names(all_datasets[["Laumer2018"]], matrix_taxa[["Laumer2018.Tplx_phylo_d1.aa"]]))$clade,
                                       extract.taxa.vector(filter.matrix.names(all_datasets[["Laumer2019"]], matrix_taxa[["Laumer2019.nonbilateria_MARE_BMGE.aa"]]))$clade,
                                       extract.taxa.vector(filter.matrix.names(all_datasets[["Moroz2014"]], matrix_taxa[["Moroz2014.ED3d.aa"]]))$clade,
                                       extract.taxa.vector(filter.matrix.names(all_datasets[["Nosenko2013"]], matrix_taxa[["Nosenko2013.nonribosomal_9187_smatrix.aa"]]))$clade,
                                       extract.taxa.vector(filter.matrix.names(all_datasets[["Nosenko2013"]], matrix_taxa[["Nosenko2013.ribosomal_14615_smatrix.aa"]]))$clade,
                                       extract.taxa.vector(all_datasets[["Philippe2009"]])$clade,
                                       extract.taxa.vector(all_datasets[["Philippe2011"]])$clade,
                                       extract.taxa.vector(all_datasets[["Pick2010"]])$clade,
                                       extract.taxa.vector(all_datasets[["Ryan2013"]])$clade,
                                       extract.taxa.vector(filter.matrix.names(all_datasets[["Simion2017"]], matrix_taxa[["Simion2017.supermatrix_97sp_401632pos_1719genes.aa"]]))$clade,
                                       extract.taxa.vector(filter.matrix.names(all_datasets[["Whelan2015"]], matrix_taxa[["Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa"]]))$clade,
                                       extract.taxa.vector(all_datasets[["Whelan2017"]])$clade),
                           "alignment" = c(rep("Borowiec0215.Best108", extract.taxa.vector(all_datasets[["Borowiec2015"]])$number),
                                           rep("Chang2015.Chang_AA", extract.taxa.vector(all_datasets[["Chang2015"]])$number),
                                           rep("Dunn2008.Dunn2008_FixedNames", extract.taxa.vector(all_datasets[["Dunn2008"]])$number),
                                           rep("Hejnol2009.Hejnol_etal_2009_FixedNames", extract.taxa.vector(all_datasets[["Hejnol2009"]])$number),
                                           rep("Laumer2018.Tplx_phylo_d1", extract.taxa.vector(filter.matrix.names(all_datasets[["Laumer2018"]], matrix_taxa[["Laumer2018.Tplx_phylo_d1.aa"]]))$number),
                                           rep("Laumer2019.nonbilateria_MARE_BMGE", extract.taxa.vector(filter.matrix.names(all_datasets[["Laumer2019"]], matrix_taxa[["Laumer2019.nonbilateria_MARE_BMGE.aa"]]))$number),
                                           rep("Moroz2014.ED3d", extract.taxa.vector(filter.matrix.names(all_datasets[["Moroz2014"]], matrix_taxa[["Moroz2014.ED3d.aa"]]))$number),
                                           rep("Nosenko2013.nonribosomal_9187_smatrix", extract.taxa.vector(filter.matrix.names(all_datasets[["Nosenko2013"]], matrix_taxa[["Nosenko2013.nonribosomal_9187_smatrix.aa"]]))$number),
                                           rep("Nosenko2013.ribosomal_14615_smatrix", extract.taxa.vector(filter.matrix.names(all_datasets[["Nosenko2013"]], matrix_taxa[["Nosenko2013.ribosomal_14615_smatrix.aa"]]))$number),
                                           rep("Philippe2009.Philippe_etal_superalignment_FixedNames", extract.taxa.vector(all_datasets[["Philippe2009"]])$number),
                                           rep("Philippe2011.UPDUNN_MB_FixedNames", extract.taxa.vector(all_datasets[["Philippe2011"]])$number),
                                           rep("Pick2010.Pick2010", extract.taxa.vector(all_datasets[["Pick2010"]])$number),
                                           rep("Ryan2013.REA_EST_includingXenoturbella", extract.taxa.vector(all_datasets[["Ryan2013"]])$number),
                                           rep("Simion2017.supermatrix_97sp_401632pos_1719genes", extract.taxa.vector(filter.matrix.names(all_datasets[["Simion2017"]], matrix_taxa[["Simion2017.supermatrix_97sp_401632pos_1719genes.aa"]]))$number),
                                           rep("Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned", extract.taxa.vector(filter.matrix.names(all_datasets[["Whelan2015"]], matrix_taxa[["Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa"]]))$number),
                                           rep("Whelan2017.Metazoa_Choano_RCFV_strict", extract.taxa.vector(all_datasets[["Whelan2017"]])$number) ) )
  # Save the data frame as a csv file
  write.csv(mastmet_df,  file = mastmet_file_path, row.names = F)
} else if (file.exists(mastmet_file_path) == TRUE){
  mastmet_df <- read.csv(mastmet_file_path, stringsAsFactors = F)
}


#### 4. Standardize the names ####
# Open the tsv files
names_df <- read.delim(taxon_table_path)
mantax_df <- read.delim(manual_taxonomy_path)
# Add column for matrix classifications for this study
names_df$MAST_matrix <- NA
# Reorganise the columns
names_df <- names_df[, c("MAST_matrix", "matrix_name", "relabelled_name", "clade_assignment", "original_matrix", "ncbi_tax_id", "ncbi_taxonomy")]
# Rename the column to match
names(names_df) <- c("MAST_matrices", "matrix_name", "relabelled_name", "clade_assignment", "Li2021_matrices", "ncbi_tax_id", "ncbi_taxonomy")
# Add in the corresponding matrix name for our study
names_df$MAST_matrix[which(names_df$Li2021_matrix == "../considered_data/Borowiec2015/Best108.nex")] <- "Borowiec2015.Best108.aa.alignment"

