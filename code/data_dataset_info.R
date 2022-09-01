## caitlinch/metazoan-mixtures/data_taxa_names.R
# Caitlin Cherryh 2022

# File to store the taxa names, clade classifications, and models of sequence evolution for the Metazoan datasets used for analysis

#### FUNCTIONS ####
count.taxa <- function(list){
  # Quick function to count the number of taxa in a list
  
  total_taxa <- c(list$Bilateria, list$Cnidaria, list$Placozoa, list$Porifera, list$Ctenophora, list$Outgroup)
  return(list("num_taxa" = length(total_taxa),
              "taxa" = total_taxa))
}

#### PHYLOGENETIC DATASETS ####
# Format for taxa names:
empty_list <- list("Bilateria" = c(),
                   "Cnidaria" = c(),
                   "Placozoa" = c(),
                   "Porifera" = c(),
                   "Ctenophora" = c(),
                   "Outgroup" = c(),
                   "Sponges_Calcarea" = c(),
                   "Sponges_Homoscleromorpha" = c(),
                   "Sponges_Hexactinellida" = c(),
                   "Sponges_Demospongiae" = c(),
                   "Sponges_1" = c(),
                   "Sponges_2" = c(),
                   "Models" = c(),
                   "Partitioned" = FALSE)

# For Dunn et. al. (2008):
dunn2008_list <- list("Bilateria" = c("Acanthoscurria_gomesiana", "Amoebidium_parasiticum", "Anoplodactylus_eroticus", "Aplysia_californica", 
                                      "Argopecten_irradians", "Asterina_pectinifera", "Biomphalaria_glabrata", "Boophilus_microplus", 
                                      "Branchiostoma_floridae", "Capitella_sp", "Carcinoscorpius_rotundicauda", "Carcinus_maenas", 
                                      "Carinoma_mutabilis", "Cerebratulus_lacteus", "Chaetoderma_nitidulum", "Chaetopleura_apiculata", 
                                      "Chaetopterus_variopedatus", "Ciona_intestinalis", "Crassostrea_virginica", "Daphnia_magna", 
                                      "Drosophila_melanogaster", "Dugesia_japonica", "Echinococcus_granulosus", "Echinoderes_horni", 
                                      "Euperipatoides_kanangrensis", "Euprymna_scolopes", "Fenneropenaeus_chinensis", "Gallus_gallus", 
                                      "Haementeria_depressa", "Homo_sapiens", "Hypsibius_dujardini", "Lumbricus_rubellus", "Macrostomum_lignano",
                                      "Mytilus_galloprovincialis", "Paraplanocera_oligoglena", "Phoronis_vancouverensis", "Platynereis_dumerilii",
                                      "Priapulus_caudatus", "Ptychodera_flava", "Richtersius_coronifer", "Saccoglossus_kowalevskii",
                                      "Schmidtea_mediterranea", "Scutigera_coleoptrata", "Spinochordodes_tellinii", "Strongylocentrotus_purpuratus",
                                      "Terebratalia_transversa", "Themiste_lageniformis", "Trichinella_spiralis", "Urechis_caupo", "Xenoturbella_bocki",
                                      "Xiphinema_index"),
                      "Cnidaria" = c("Acropora_millepora", "Nematostella_vectensis", "Cyanea_capillata", "Hydra_magnipapillata", "Hydractinia_echinata"),
                      "Placozoa" = c(),
                      "Porifera" = c("Oscarella_carmela"),
                      "Ctenophora" = c("Mertensiid_sp", "Mnemiopsis_leidyi"),
                      "Outgroup" = c("Monosiga_ovata", "Capsaspora_owczarzaki", "Sphaeroforma_arctica", "Saccharomyces_cerevisiae", "Cryptococcus_neoformans"),
                      "Sponges_Calcarea" = c(),
                      "Sponges_Homoscleromorpha" = c("Oscarella_carmela"),
                      "Sponges_Hexactinellida" = c(),
                      "Sponges_Demospongiae" = c(),
                      "Sponges_1" = c("Sponges_Calcarea", "Sponges_Homoscleromorpha"),
                      "Sponges_2" = c("Sponges_Hexactinellida", "Sponges_Demospongiae"),
                      "Models" = c(),
                      "Partitioned" = FALSE)

# For Hejnol et. al. (2009):
hejnol2009_list <- list("Bilateria" = c(),
                        "Cnidaria" = c(),
                        "Placozoa" = c(),
                        "Porifera" = c(),
                        "Ctenophora" = c(),
                        "Outgroup" = c(),
                        "Sponges_Calcarea" = c(),
                        "Sponges_Homoscleromorpha" = c(),
                        "Sponges_Hexactinellida" = c(),
                        "Sponges_Demospongiae" = c(),
                        "Sponges_1" = c(),
                        "Sponges_2" = c(),
                        "Models" = c(),
                        "Partitioned" = FALSE)

# For Whelan et. al. (2017):
whelan2017_list <- list("Bilateria" = c("Homo_sapiens", "Strongylocentrotus_purpatus", "Hemithris_psittacea", "Capitella_teleta", "Drosophila_melanogaster","Daphnia_pulex"),
                        "Cnidaria" = c("Hydra_vulgaris", "Bolocera_tuediae", "Aiptasia_pallida", "Hormathia_digitata", "Nematostella_vectensis", "Acropora_digitifera", 
                                       "Eunicella_verrucosa", "Hydra_viridissima", "Hydra_oligactis", "Physalia_physalia", "Abylopsis_tetragona","Craseo_lathetica",
                                       "Nanomia_bijuga", "Agalma_elegans", "Periphyla_periphyla"),
                        "Placozoa" = c("Trichoplax_adhaerens"),
                        "Porifera" = c("Cliona_varians", "Sycon_coactum", "Sycon_ciliatum", "Corticium_candelabrum", "Oscarella_carmela", "Hyalonema_populiferum",
                                       "Aphrocallistes_vastus", "Rossella_fibulata", "Sympagella_nux", "Ircinia_fasciculata", "Chondrilla_nucula", "Amphimedon_queenslandica",
                                       "Petrosia_ficiformis", "Spongilla_lacustris", "Pseudospongosorites_suberitoides", "Mycale_phylophylla", "Latrunculia_apicalis", 
                                       "Crella_elegans", "Kirkpatrickia_variolosa"),
                        "Ctenophora" = c("Euplokamis_dunlapae", "Vallicula_sp", "Coeloplana_astericola", "Hormiphora_californica", "Hormiphora_palmata", "Pleurobrachia_pileus",
                                         "Pleurobrachia_bachei", "Pleurobrachia_sp_South_Carolina_USA", "Cydippida_sp_Maryland_USA", "Callianira_Antarctica", "Mertensiidae_sp_Antarctica",
                                         "Mertensiidae_sp_Washington_USA", "Cydippida_sp", "Dryodora_glandiformis", "Lobatolampea_tetragona", "Beroe_abyssicola", "Beroe_sp_Antarctica",
                                         "Beroe_ovata", "Beroe_sp_Queensland_Australia", "Beroe_forskalii", "Ocyropsis_sp_Bimini_Bahamas", "Ocyropsis_crystallina", "Ocyropsis_sp_Florida_USA",
                                         "Bolinopsis_infundibulum", "Mnemiopsis_leidyi", "Bolinopsis_ashleyi", "Lobata_sp_Punta_Arenas_Argentina", "Eurhamphaea_vexilligera", "Cestum_veneris",
                                         "Ctenophora_sp_Florida_USA"),
                        "Outgroup" = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"),
                        "Sponges_Calcarea" = c("Sycon_coactum", "Sycon_ciliatum"),
                        "Sponges_Homoscleromorpha" = c("Oscarella_carmela", "Corticium_candelabrum"),
                        "Sponges_Hexactinellida" = c("Hyalonema_populiferum", "Sympagella_nux", "Rossella_fibulata", "Aphrocallistes_vastus"),
                        "Sponges_Demospongiae" = c("Ircinia_fasciculata", "Chondrilla_nucula", "Spongilla_lacustris", "Cliona_varians", "Pseudospongosorites_suberitoides", 
                                                   "Mycale_phylophylla", "Latrunculia_apicalis", "Kirkpatrickia_variolosa", "Crella_elegans", "Petrosia_ficiformis", "Amphimedon_queenslandica"),
                        "Sponges_1" = c("Sponges_Calcarea", "Sponges_Homoscleromorpha"),
                        "Sponges_2" = c("Sponges_Hexactinellida", "Sponges_Demospongiae"),
                        "Models" = c("PartitionFinder"),
                        "Partitioned" = TRUE)

# Create one list that contains all other lists, indexed by dataset identifier (first author+year of publication)
all_datasets <- list("Dunn2008" = dunn2008_list,
                     "Hejnol2009" = hejnol2009_list,
                     "Whelan2017" = whelan2017_list)



#### MODELS ####
redmond_models <- c("C30", "C40", "C50", "C60", "EHO+G", "EX_EHO+G", "EX2+G", "EX3+G", 
                    "JTT+C20", "JTT+C30", "JTT+C40", "JTT+C50", "JTT+C60", "JTT+G4", "JTTDCMut+G4", "JTTDCMut+I+G4", 
                    "LG+C30", "LG+C40", "LG+C50", "LG+C60", "LG+CF4", "LG+F+G4", "LG+F+I+G4", "LG+G4", "LG+I+G4", "LG4M", 
                    "mtZOA+F+G4", "mtZOA+G4", "mtZOA+I+G4", "PMB+G4", "rtREV+G4", "rtREV+I+G4", "UL2+G", "UL3+G", 
                    "WAG+C20", "WAG+C30", "WAG+C40", "WAG+C50", "WAG+C60", "WAG+G4", "WAG+I+G4")
li_models <- c("C60+LG", "C60+Poisson", "C60+WAG", "CAT+F81", "CAT+GTR", "CAT+GTR+G", "CAT+Poisson",
               "CAT+Poisson+G", "CAT+WAG+F", "GTR+FO", "GTR+G+FO", "GTR20", "LG+G+F", "WAG",        
               "WAG+F", "WAG+G+F" )
spreadsheet_models <- c("ModelFinder", "rtREV+G+F", "GTR+G4", 
                        "WAG", "WAG+F", "WAG+CAT+G4", "WAG+C60", "WAG+C50", "WAG+C40", "WAG+C30", "WAG+C20", "WAG+C10",
                        "LG", "LG+F", "LG+G", "LG+G4+F", "LG+C60", "LG+C50", "LG+C40", "LG+C30", "LG+C20", "LG+C10",
                        "C20+LG+FO+R4", "C60+LG+FO+R4", "C60+LG+G+F", "LG+PMSF+G", 
                        "Poisson+G", "Poisson+CAT", "Poisson+CAT+G", "JTT", "JTT+C60", "JTT+C50", "JTT+C40", "JTT+C30", "JTT+C20", "JTT+C10",
                        "EX2+G", "EX3+G", "EX_EHO+G", "C60", "C50", "C40", "C30", "C20", "C10", 
                        "UL3+G", "UL2+G", "mtZOA+G4", "LG4M")
# Collate all models
all_models <- c(redmond_models, li_models, spreadsheet_models)
# Sort and identify unique models
all_models <- sort(unique(all_models))


