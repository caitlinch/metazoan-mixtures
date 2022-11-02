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
convert.alignment.name <- function(dataset_name, alignment_name){
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

check.manual.taxonomy.map <- function(species_row, manual_taxonomy_df){
  # Function to take a species name, and check if it's in the manual taxonomy map
  
  # Remove any strings of underscores or trailing underscores from the species name
  reformat_species_name <- gsub("_$","",gsub("\\_\\_", "", species_row$original_name))
  # Check for identical names
  matching_ind <- grep(reformat_species_name, manual_taxonomy_df$original_name)
  # Check whether there is a matching name 
  if (identical(matching_ind, integer(0)) == FALSE){
    matching_name <- manual_taxonomy_df$new_name[matching_ind]
  } else if (identical(matching_ind, integer(0)) == TRUE){
    # There is no matching name
    # Return NA
    matching_name <- NA
  }
  
  # Return the matched name
  return(matching_name)
}


find.species.name <- function(species_row, taxon_table_df, manual_taxonomy_df){
  # Function to take a species name and check Li et. al. tsv files to see if a reconciled species name exists
  # Works for all datasets except Laumer2018, Laumer2019, Philippe2011, Pick2010, Simion2017
  
  # Determine which Li et. al. 2021 alignment corresponds to this alignment
  li_alignment_name <- convert.alignment.name(species_row$dataset, species_row$alignment)
  # Check whether there is a corresponding Li et. al. alignment
  if (is.na(li_alignment_name) == FALSE){
    # If there is a corresponding alignment in Li et. al., find the relabelled species name for this species
    # Reduce the taxon_table_df to just the species present in this dataset
    species_df <- taxon_table_df[(intersect(grep(species_row$dataset, taxon_table_df$original_matrix), grep(li_alignment_name, taxon_table_df$original_matrix))),]
    # Process the species data frame, if required
    if (species_row$dataset == "Philippe2009" & li_alignment_name == "Philippe2009"){
      # If this is the Philippe2009 dataset, I manually fixed the taxa names to remove the underscores ("____")
      # Remove the strings of underscores from the end of the species names
      species_df$matrix_name <- gsub("\\_\\_", "", species_df$matrix_name)
      # If there is a single trailing underscore left, remove it
      species_df$matrix_name <- gsub("_$","",species_df$matrix_name)
    }
    # Check whether this species name is present in the tsv files
    relabelled_name <- species_df$relabelled_name[which(species_df$matrix_name == species_row$original_name)]
  } else if (is.na(li_alignment_name) == TRUE){
    # There is no corresponding alignment
    # This species will need further checks
    if (species_row$dataset == "Philippe2011"){
      # For Philippe2011 dataset
      # Check if any of the matrix names have this name
      check_df <- taxon_table_df[taxon_table_df$matrix_name == species_row$original_name, ]
      # Check that all taxa with this matrix name were given the same original name
      if (length(unique(check_df$relabelled_name)) == 1){
        # If the relabelled names from this matrix name are identical, return the relabelled name
        relabelled_name <- unique(check_df$relabelled_name)
      } else {
        # Check the list of manual conversions to see if there's a matching name
        relabelled_name <- check.manual.taxonomy.map(species_row, manual_taxonomy_df)
      }
    } else {
      # Return NA
      relabelled_name <- NA 
    } # end if (species_row$dataset == "Philippe2011")
  } # end if (is.na(li_alignment_name) == FALSE)
  
  # Return the relabelled species name
  return(relabelled_name)
}








