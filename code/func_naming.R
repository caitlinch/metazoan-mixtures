## caitlinch/metazoan-mixtures/func_names.R
# Caitlin Cherryh 2022

# Functions for taxa naming

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




