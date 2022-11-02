## caitlinch/metazoan-mixtures/func_names.R
# Caitlin Cherryh 2022

# Functions for taxa naming


#### Functions for collecting and summarising information from the dataset taxa lists ### 
extract.taxa.vector <- function(dataset_list){
  # Given a taxa list, this function returns all taxa in that dataset and the number of taxa in that dataset
  
  # Collate all the taxa
  if ("Placozoa" %in% names(dataset_list) == TRUE){
    dl_taxa <- c(dataset_list$Bilateria, dataset_list$Cnidaria, dataset_list$Placozoa, 
                 dataset_list$Porifera, dataset_list$Ctenophora, dataset_list$Outgroup)
    dl_classifications <- c(rep("Bilateria", length(dataset_list$Bilateria)), rep("Cnidaria", length(dataset_list$Cnidaria)), 
                            rep("Placozoa", length(dataset_list$Placozoa)), rep("Porifera", length(dataset_list$Porifera)),
                            rep("Ctenophora", length(dataset_list$Ctenophora)), rep("Outgroup", length(dataset_list$Outgroup)))
    dl_num_taxa <- length(dl_taxa)
  } else if ("Placozoa" %in% names(dataset_list) == FALSE){
    dl_taxa <- c(dataset_list$Bilateria, dataset_list$Cnidaria, dataset_list$Porifera, 
                 dataset_list$Ctenophora, dataset_list$Outgroup)
    dl_classifications <- c(rep("Bilateria", length(dataset_list$Bilateria)), rep("Cnidaria", length(dataset_list$Cnidaria)), 
                            rep("Porifera", length(dataset_list$Porifera)), rep("Ctenophora", length(dataset_list$Ctenophora)), 
                            rep("Outgroup", length(dataset_list$Outgroup)))
    dl_num_taxa <- length(dl_taxa)
  } # end if if ("Placozoa" %in% names(dataset_list) == TRUE)
  
  # Create the output list
  output_list = list(taxa = dl_taxa, clade = dl_classifications, number = dl_num_taxa)
  
  # Return output
  return(output_list)
}


filter.matrix.names <- function(taxa_list, matrix_subset){
  # Function to remove all names from a taxa list that are not present in the matrix of interest
  # Returns a taxa list (with only taxa in the alignments present)
  
  # Create a new taxa list (removing taxa as you go)
  new_taxa_list <- list("Bilateria" = taxa_list$Bilateria[taxa_list$Bilateria %in% matrix_subset],
                        "Cnidaria" = taxa_list$Cnidaria[taxa_list$Cnidaria %in% matrix_subset],
                        "Placozoa" = taxa_list$Placozoa[taxa_list$Placozoa %in% matrix_subset],
                        "Porifera" = taxa_list$Porifera[taxa_list$Porifera %in% matrix_subset],
                        "Ctenophora" = taxa_list$Ctenophora[taxa_list$Ctenophora %in% matrix_subset],
                        "Outgroup" = taxa_list$Outgroup[taxa_list$Outgroup %in% matrix_subset],
                        "Sponges_Calcarea" = taxa_list$Sponges_Calcarea[taxa_list$Sponges_Calcarea %in% matrix_subset],
                        "Sponges_Homoscleromorpha" = taxa_list$Sponges_Homoscleromorpha[taxa_list$Sponges_Homoscleromorpha %in% matrix_subset],
                        "Sponges_Hexactinellida" = taxa_list$Sponges_Hexactinellida[taxa_list$Sponges_Hexactinellida %in% matrix_subset],
                        "Sponges_Demospongiae" = taxa_list$Sponges_Demospongiae[taxa_list$Sponges_Demospongiae %in% matrix_subset],
                        "Sponges_1" = taxa_list$Sponges_1,
                        "Sponges_2" = taxa_list$Sponges_2,
                        "Models" = taxa_list$Models,
                        "Partitioned" = taxa_list$Partitioned,
                        "Estimate.Paraphyletic.Sponges" = taxa_list$Estimate.Paraphyletic.Sponges
                        )
  # Output the filtered list
  return(new_taxa_list)
}


#### Functions to reconcile species names across datasets ####
match.alignment.name <- function(dataset_name, alignment_name){
  # Small function to take an alignment and dataset and return the corresponding alignment name 
  #   for the taxon_table_df from Li et. al. (2021)
  
  if (dataset_name == "Borowiec2015" & alignment_name == "Best108"){
    li_name <- "Best108"
  } else if (dataset_name == "Chang2015" & alignment_name == "Chang_AA"){
    li_name <- "Chang2015"
  } else if (dataset_name == "Dunn2008" & alignment_name == "Dunn2008_FixedNames"){
    li_name <- "Dunn2008"
  } else if (dataset_name == "Hejnol2009" & alignment_name == "Hejnol_etal_2009_FixedNames"){
    li_name <- "Hejnol2009"
  } else if (dataset_name == "Moroz2014" & alignment_name == "ED3d"){
    li_name <- "3d"
  } else if (dataset_name == "Nosenko2013" & alignment_name == "nonribosomal_9187_smatrix"){
    li_name <- "nonribo_9187"
  } else if (dataset_name == "Nosenko2013" & alignment_name == "ribosomal_14615_smatrix"){
    li_name <- "ribo_11057"
  } else if (dataset_name == "Philippe2009" & alignment_name == "Philippe_etal_superalignment_FixedNames"){
    li_name <- "Philippe2009"
  } else if (dataset_name == "Ryan2013" & alignment_name == "REA_EST_includingXenoturbella"){
    li_name <- "est"
  } else if (dataset_name == "Whelan2015" & alignment_name == "Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned"){
    li_name <- "D10"
  } else if (dataset_name == "Whelan2017" & alignment_name == "Metazoa_Choano_RCFV_strict"){
    li_name <- "full"
  } else {
    li_name <- NA
  }
  
  # Return the corresponding alignment name from Li et. al. (2021)
  return(li_name)
}

find.species.name <- function(species_row, taxon_table_df, manual_taxonomy_df){
  # Function to take a species name and check Li et. al. tsv files to see if a reconciled species name exists
  
  # Determine which Li et. al. 2021 alignment corresponds to this alignment
  li_alignment_name <- match.alignment.name(species_row$dataset, species_row$alignment)
  # Reduce the taxon_table_df to just the species present in this dataset
  
  # Check whether this species name is present in the tsv files
  
}








