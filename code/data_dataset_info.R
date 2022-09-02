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

check.sponges <- function(list){
  # Quick function to make sure that all the sponge species are included in the correct places in the list
  
  # Extract species
  Porifera_taxa <- list$Porifera
  combined_taxa <- c(list$Sponges_Calcarea, list$Sponges_Homoscleromorpha, 
                     list$Sponges_Hexactinellida, list$Sponges_Demospongiae)
  # Check for equality
  sponge_check <- setequal(Porifera_taxa, combined_taxa)
  # Collate output
  op_list <- list("Matching" = sponge_check,
                  "Porifera_length" = length(Porifera_taxa),
                  "Porifera_taxa" = Porifera_taxa,
                  "Combined_length" = length(combined_taxa),
                  "Combined_taxa" = combined_taxa)
  return(op_list)
}

check.outgroups <- function(list){
  # Quick function to make sure all the outgroup species are included in the correct places on the list
  
  # Extract species
  all_outgroup_taxa <- list$Outgroup
  collated_outgroup_taxa <- unique(unlist(list[grep("Outgroup_", names(list), value = T)]))
  # Check for equality
  outgroup_check <- setequal(all_outgroup_taxa, collated_outgroup_taxa)
  # Collate output
  op_list <- list("Matching" = outgroup_check,
                  "Outgroup_length" = length(all_outgroup_taxa),
                  "Outgroup_taxa" = all_outgroup_taxa,
                  "Combined_length" = length(collated_outgroup_taxa),
                  "Combined_taxa" = collated_outgroup_taxa)
  return(op_list)
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
# Note: Taxa "Convoluta_pulchra" (a marine worm) is labelled "Isodiametra_pulchra" on Figure 2, the main phylogeny for this paper. The current accepted name 
#       for this taxa is "Aphanostoma pulchra".
hejnol2009_list <- list("Bilateria" = c("Xenoturbella_bocki", "Meara_stichopi", "Nemertoderma_westbladi", "Neochildia_fusca", "Symsagittifera_roscoffensis",
                                        "Convolutriloba_longifissura", "Convoluta_pulchra", "Priapulus_caudatus", "Spinochordodes_tellinii", "Echinoderes_horni",
                                        "Xiphinema_index", "Trichinella_spiralis", "Hypsibius_dujardini", "Richtersius_coronifer", "Epiperipatus_sp", 
                                        "Euperipatoides_kanangrensis", "Drosophila_melanogaster", "Onychiurus_arcticus", "Daphnia_pulex", "Carcinus_maenas",
                                        "Fenneropenaeus_chinensis", "Scutigera_coleoptrata", "Anoplodactylus_eroticus", "Acanthoscurria_gomesiana", 
                                        "Boophilus_microplus", "Carcinoscorpius_rotundicauda", "Phoronis_vancouverensis", "Spadella_cephaloptera",
                                        "Flaccisagitta_enflata", "Platynereis_dumerilii", "Capitella_sp", "Urechis_caupo", "Myzostoma_seymourcollegiorum",
                                        "Lumbricus_rubellus", "Haementeria_depressa", "Helobdella_robusta", "Themiste_lageniformis", "Chaetopterus_sp",
                                        "Chaetopleura_apiculata", "Euprymna_scolopes", "Chaetoderma_nitidulum", "Lottia_gigantea", "Mytilus_galloprovincialis",
                                        "Crassostrea_virginica", "Argopecten_irradians", "Aplysia_californica", "Biomphalaria_glabrata", "Terebratalia_transversa",
                                        "Cerebratulus_lacteus", "Carinoma_mutabilis", "Gnathostomula_peregrina", "Turbanella_ambronensis", "Macrostomum_lignano",
                                        "Paraplanocera_sp", "Taenia_solium", "Echinococcus_granulosus", "Dugesia_japonica", "Schmidtea_mediterranea", "Philodina_roseola",
                                        "Brachionus_plicatilis", "Bugula_neritina", "Cristatella_mucedo", "Symbion_pandora", "Pedicellina_sp", "Pedicellina_cernua",
                                        "Branchiostoma_floridae", "Ciona_intestinalis", "Diplosoma_listerianum", "Halocynthia_roretzi", "Gallus_gallus", "Homo_sapiens",
                                        "Saccoglossus_kowalevskii", "Ptychodera_flava", "Strongylocentrotus_purpuratus", "Asterina_pectinifera"),
                        "Cnidaria" = c("Acropora_millepora", "Nematostella_vectensis", "Cyanea_capillata", "Hydra_magnipapillata", "Hydractinia_echinata"),
                        "Placozoa" = c("Trichoplax_adhaerens"),
                        "Porifera" = c("Amphimedon_queenslandica", "Suberites_domuncula", "Oscarella_carmela"),
                        "Ctenophora" = c("Mertensiid_sp", "Mnemiopsis_leidyi", "Pleurobrachia_pileus"),
                        "Outgroup" = c("Saccharomyces_cerevisiae", "Cryptococcus_neoformans", "Capsaspora_owczarzaki", "Sphaeroforma_arctica", "Amoebidium_parasiticum",
                                       "Monosiga_ovata",  "Monosiga_brevicollis"),
                        "Sponges_Calcarea" = c(),
                        "Sponges_Homoscleromorpha" = c("Oscarella_carmela"),
                        "Sponges_Hexactinellida" = c(),
                        "Sponges_Demospongiae" = c("Suberites_domuncula", "Amphimedon_queenslandica"),
                        "Sponges_1" = c("Sponges_Calcarea", "Sponges_Homoscleromorpha"),
                        "Sponges_2" = c("Sponges_Hexactinellida", "Sponges_Demospongiae"),
                        "Models" = c(),
                        "Partitioned" = FALSE)

# For Philippe et. al. (2009):
# Note: Taxa "RenieraXsp" (a marine sponge) is labelled "Amphimedon" in Figure 1, the main phylogeny for this paper. The current accepted name for this 
#       species is "Amphimedon queenslandica".
# Note: One species not shown in main figure (Figure 1) but in alignment: Euperipatoides kanangrensis (a velvet worm from Bilateria)
philippe2009_list <- list("Bilateria" = c("Xenoturbella_bocki", "Strongylocentrotus_purpuratus", "Saccoglossus_kowalevskii", "Branchiostoma_floridae", "Petromyzon_marinus",
                                          "Danio_rerio", "Molgula_tectiformis", "CionaXsavignyi", "Euprymna_scolopes", "Crassostrea_gigas", "AplysiaXcalifornica", 
                                          "Pedicellina_cernua", "Capitella_sp__i_ecs-2004", "TubifexXtubifex", "Helobdella_robusta", "DaphniaXpulex", "PediculusXhumanus",
                                          "NasoniaXvitripennis", "ScutigeraXcoleoptrata", "IxodesXscapularis", "Anoplodactylus_eroticus"),
                          "Cnidaria" = c("Nematostella_vectensis", "Metridium_senile", "Montastraea_faveolata", "Acropora_millepora", "CyaneaXcapillata", "HydraXmagnipapillata",
                                         "ClytiaXhemisphaerica", "Podocoryne_carnea", "Hydractinia_echinata", "Euperipatoides_kanangrensis"),
                          "Placozoa" = c("Trichoplax_adhaerens"),
                          "Porifera" = c("Oscarella_carmela", "SyconXraphanus", "LeucettaXchagosensis", "OopsacasXminuta", "Heterochone_sp", "Carteriospongia_foliascens",
                                         "RenieraXsp", "Suberites_domuncula", "Ephydatia_muelleri"),
                          "Ctenophora" = c("Pleurobrachia_pileus", "mertensiid_sp", "Mnemiopsis_leidyi"),
                          "Outgroup" = c("Spizellomyces_punctatus", "Batrachochytrium_dendrobatidis", "AllomycesXmacrogynus", "RhizopusXoryzae", "Phycomyces_blakesleeanus",
                                         "Sphaeroforma_arctica", "Amoebidium_parasiticum", "Capsaspora_owczarzaki", "MonosigaXovata", "Proterospongia_sp", "MonosigaXbrevicollis"),
                          "Outgroup_Choanoflagellata" = c("MonosigaXovata", "Proterospongia_sp", "MonosigaXbrevicollis"),
                          "Outgroup_Fungi" = c("Spizellomyces_punctatus", "Batrachochytrium_dendrobatidis", "AllomycesXmacrogynus", "RhizopusXoryzae", "Phycomyces_blakesleeanus",
                                               "Sphaeroforma_arctica", "Amoebidium_parasiticum", "Capsaspora_owczarzaki"),
                          "Sponges_Calcarea" = c("SyconXraphanus", "LeucettaXchagosensis"),
                          "Sponges_Homoscleromorpha" = c("Oscarella_carmela"),
                          "Sponges_Hexactinellida" = c("OopsacasXminuta", "Heterochone_sp"),
                          "Sponges_Demospongiae" = c("Carteriospongia_foliascens", "RenieraXsp", "Suberites_domuncula", "Ephydatia_muelleri"),
                          "Sponges_1" = c("Sponges_Calcarea", "Sponges_Homoscleromorpha"),
                          "Sponges_2" = c("Sponges_Hexactinellida", "Sponges_Demospongiae"),
                          "Models" = c(),
                          "Partitioned" = FALSE)

# For Pick et. al. (2010):
pick2010_list <- list("Bilateria" = c("Lumbricus", "Haementeri", "Urechis_ca", "Capitella", "Platynerei", "Themiste_l", "Chaetopter", "Terebratal", "Phoronis_v", "Cerebratul",
                                      "Carinoma_m", "Crassostre", "Argopecten", "Mytilus_ga", "Biomphalar", "Aplysia_ca", "Euprymna_s", "Chaetopleu", "Chaetoderm", "Schmidtea",
                                      "Dugesia_ja", "Echinococc", "Paraplanoc", "Macrostomu", "Acanthoscu", "Boophilus", "Anoplodact", "Carcinosco", "Scutigera", "Litopenaeu",
                                      "Homarus_am", "Drosophila", "Daphnia_pu", "Euperipato", "Xiphinema", "Trichinell", "Spinochord", "Richtersiu", "Hypsibius", "Priapulus",
                                      "Echinodere", "Saccogloss", "Ptychodera", "Strongyloc", "Asterina_p", "Xenoturbel", "Homo_sapie", "Gallus_gal", "Ciona_inte", "Branchiost"),
                      "Cnidaria" = c("Podocoryne", "Hydractini", "Hydra_magn", "Clytia_hem", "Cyanea_cap", "Metridium", "Anemonia_v", "Nematostel", "Montastrae", "Acropora_m"),
                      "Placozoa" = c("Trichoplax"),
                      "Porifera" = c("Lubomirski", "Ephydatia", "Pachydicty", "Suberites", "Amphimedon", "Carteriosp", "Heterochon", "Crateromor", "Oopsacas_m", "Oscarella",
                                     "Oscarell00", "Sycon_raph", "Leucetta_c"),
                      "Ctenophora" = c("mertensiid", "Pleurobrac", "Mnemiopsis"),
                      "Outgroup" = c("Proterospo", "Monosiga_b", "Monosiga00", "Amoebidium", "Capsaspora", "Sphaerofor"),
                      "Outgroup_Choanoflagellata" = c("Proterospo", "Monosiga_b", "Monosiga00"),
                      "Outgroup_Holozoa" = c("Amoebidium", "Capsaspora", "Sphaerofor"),
                      "Outgroup_Opisthokont" = c("Amoebidium", "Capsaspora", "Sphaerofor"),
                      "Sponges_Calcarea" = c("Sycon_raph", "Leucetta_c"),
                      "Sponges_Homoscleromorpha" = c("Oscarella", "Oscarell00"),
                      "Sponges_Hexactinellida" = c("Heterochon", "Crateromor", "Oopsacas_m"),
                      "Sponges_Demospongiae" = c("Lubomirski", "Ephydatia", "Pachydicty", "Suberites", "Amphimedon", "Carteriosp"),
                      "Sponges_1" = c("Sponges_Calcarea", "Sponges_Homoscleromorpha"),
                      "Sponges_2" = c("Sponges_Hexactinellida", "Sponges_Demospongiae"),
                      "Models" = c(),
                      "Partitioned" = FALSE)

# For Philippe et. al. (2011):
# Note: Taxa "Reniera_sp" (a marine sponge) is labelled "Amphimedon" in Figure 1, the main phylogeny for this paper. The current accepted name for this 
#       species is "Amphimedon queenslandica".
philippe2011_list <- list("Bilateria" = c("Acanthoscurria_gomesiana", "Amoebidium_parasiticum", "Anoplodactylus_eroticus", "Aplysia_californica", 
                                          "Argopecten_irradians", "Asterina_pectinifera", "Biomphalaria_glabrata", "Boophilus_microplus", 
                                          "Branchiostoma_floridae", "Carcinoscorpius_rotundicauda", "Carcinus_maenas", "Carinoma_mutabilis",
                                          "Cerebratulus_lacteus", "Chaetoderma_nitidulum", "Chaetopleura_apiculata", "Ciona_intestinalis", 
                                          "Drosophila_melanogaster", "Dugesia_japonica", "Echinococcus_granulosus", "Echinoderes_horni", 
                                          "Euperipatoides_kanangrensis", "Euprymna_scolopes", "Fenneropenaeus_chinensis", "Gallus_gallus", 
                                          "Haementeria_depressa", "Homo_sapiens", "Hypsibius_dujardini", "Lumbricus_rubellus", "Macrostomum_lignano",
                                          "Mytilus_galloprovincialis", "Phoronis_vancouverensis", "Platynereis_dumerilii", "Priapulus_caudatus",
                                          "Ptychodera_flava", "Richtersius_coronifer", "Saccoglossus_kowalevskii", "Schmidtea_mediterranea", 
                                          "Scutigera_coleoptrata", "Spinochordodes_tellinii", "Strongylocentrotus_purpuratus", "Terebratalia_transversa",
                                          "Themiste_lageniformis", "Trichinella_spiralis", "Urechis_caupo", "Xenoturbella_bocki", "Xiphinema_index",
                                          "Daphnia_pulex", "Crassostrea_gigas", "Capitella_sp__i_ecs-2004", "Chaetopterus_sp", "Pedicellina_cernua",
                                          "Paraplanocera_sp", "Turbanella_ambronensis", "Myzostoma_seymourcollegiorum", "Symsagittifera_roscoffensis",
                                          "Neochildia_fusca", "Gnathostomula_peregrina", "Philodina_roseola", "Brachionus_plicatilis", "Bugula_neritina",
                                          "Cristatella_mucedo", "Spadella_cephaloptera", "Flaccisagitta_enflata"),
                          "Cnidaria" = c("Hydractinia_echinata", "Hydra_magnipapillata", "Cyanea_capillata", "Nematostella_vectensis", "Acropora_millepora"),
                          "Placozoa" = c(),
                          "Porifera" = c("Oscarella_carmela", "Reniera_sp"),
                          "Ctenophora" = c("Mnemiopsis_leidyi", "mertensiid_sp"),
                          "Outgroup" = c("Saccharomyces_cerevisiae", "Cryptococcus_neoformans", "Sphaeroforma_arctica", "Amoebidium_parasiticum", 
                                         "Capsaspora_owczarzaki", "Monosiga_ovata"),
                          "Sponges_Calcarea" = c(),
                          "Sponges_Homoscleromorpha" = c("Oscarella_carmela"),
                          "Sponges_Hexactinellida" = c(),
                          "Sponges_Demospongiae" = c("Reniera_sp"),
                          "Sponges_1" = c("Sponges_Calcarea", "Sponges_Homoscleromorpha"),
                          "Sponges_2" = c("Sponges_Hexactinellida", "Sponges_Demospongiae"),
                          "Models" = c(),
                          "Partitioned" = FALSE)

# For Nosenko et. al. (2013):
nosenko2013_list <- list("Bilateria" = c("Petromy_ma", "Danio_reri", "Molgula_te", "Ciona_inte", "Branchi_fl", "Saccogl_ko", "Tubifex_tu", "Helobde_ro", "Capitel_sp",
                                         "Crassos_gi", "Aplysia_ca", "Nasonia_vi", "Daphnia_pu", "Ixodes_sca", "Anoplod_er", "Euprymn_sc", "Pedicel_ce", "Pedicul_hu",
                                         "Scutige_co", "Strongy_pu", "Xenotur_bo"),
                         "Cnidaria" = c("Podocor_ca", "Hydract_ec", "Hydra_vulg", "Hydra_magn", "Clytia_hem", "Nematos_ve", "Metridi_se", "Acropor_mi", "Acropor_hy",
                                        "Acropor_pa", "Porites", "Cyanea_cap", "Montast_fa"),
                         "Placozoa" = c("Trichop_ad", "Trichop_H4"),
                         "Porifera" = c("Sycon_cili", "Sycon_raph", "Leucett_ch", "Leucoso_sp", "Oscarel_sp", "Oscarel_lo", "Oscarel_ca", "Corticium", "Cratero_me",
                                        "Oopsacas_m", "Heteroc_ca", "Aphroca_va", "Suberit_do", "Suberit_fu", "Asbesto_sp", "Tethya_wil", "Astrosc_wi", "Lubomir_ba",
                                        "Ephydat_mu", "Pachydi_gl", "Amphime_qu", "Carteri_fo"),
                         "Ctenophora" = c("Mnemiop_le", "Beroe_sp", "Mertens_sp", "Pleurob_pi"),
                         "Outgroup" = c("Monosig_ov", "Protero_sp", "Monosig_br", "Allomyc_ma", "Batrach_de", "Phycomy_bl", "Spizell_pu", "Capsasp_ow", "Amoebid_pa",
                                        "Sphaero_ar"),
                         "Outgroup_Choanoflagellata" = c("Monosig_ov", "Protero_sp", "Monosig_br"),
                         "Outgroup_Fungi" = c("Allomyc_ma", "Batrach_de", "Phycomy_bl", "Spizell_pu"), 
                         "Outgroup_Filasterea" = c("Capsasp_ow"),
                         "Outgroup_Ichthyosporea" = c("Amoebid_pa", "Sphaero_ar"), 
                         "Sponges_Calcarea" = c("Sycon_cili", "Sycon_raph", "Leucett_ch", "Leucoso_sp"),
                         "Sponges_Homoscleromorpha" = c("Oscarel_sp", "Oscarel_lo", "Oscarel_ca", "Corticium"),
                         "Sponges_Hexactinellida" = c("Cratero_me", "Oopsacas_m", "Heteroc_ca", "Aphroca_va"),
                         "Sponges_Demospongiae" = c("Suberit_do", "Suberit_fu", "Asbesto_sp", "Tethya_wil", "Astrosc_wi", "Lubomir_ba", "Ephydat_mu", "Pachydi_gl",
                                                    "Amphime_qu", "Carteri_fo"),
                         "Sponges_1" = c("Sponges_Calcarea", "Sponges_Homoscleromorpha"),
                         "Sponges_2" = c("Sponges_Hexactinellida", "Sponges_Demospongiae"),
                         "Models" = c(),
                         "Partitioned" = FALSE)

# For Ryan et. al. (2013):
ryan2013_list <- list("Bilateria" = c(),
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
                     "Philippe2009" = philippe2009_list,
                     "Pick2010" = pick2010_list,
                     "Philippe2011" = philippe2011_list,
                     "Nosenko2013" = nosenko2013_list,
                     "Ryan2013" = ryan2013_list,
                     "Whelan2017" = whelan2017_list)

# Create a list that contains the taxa within each matrix
#   For use when there are multiple alignments for a single dataset identification
#   Checks this list to get the set of taxa to use for estimating constraint trees
matrix_taxa <- list("Nosenko2013.nonribosomal_9187_smatrix.aa" = c("Acropor_hy", "Acropor_mi", "Acropor_pa", "Allomyc_ma", "Amoebid_pa", "Amphime_qu", "Anoplod_er",
                                                                   "Aplysia_ca", "Aphroca_va", "Asbesto_sp", "Astrosc_wi", "Batrach_de", "Beroe_sp", "Branchi_fl",
                                                                   "Capitel_sp", "Capsasp_ow", "Carteri_fo", "Ciona_inte", "Clytia_hem", "Corticium", "Crassos_gi",
                                                                   "Cratero_me", "Cyanea_cap", "Danio_reri", "Daphnia_pu", "Ephydat_mu", "Euprymn_sc", "Helobde_ro",
                                                                   "Heteroc_ca", "Hydra_magn", "Hydra_vulg", "Hydract_ec", "Ixodes_sca", "Leucett_ch", "Leucoso_sp",
                                                                   "Lubomir_ba", "Mertens_sp", "Metridi_se", "Mnemiop_le", "Molgula_te", "Monosig_br", "Monosig_ov",
                                                                   "Montast_fa", "Nasonia_vi", "Nematos_ve", "Oopsacas_m", "Oscarel_ca", "Oscarel_lo", "Oscarel_sp",
                                                                   "Pachydi_gl", "Pedicel_ce", "Pedicul_hu", "Petromy_ma", "Phycomy_bl", "Pleurob_pi", "Podocor_ca",
                                                                   "Porites", "Protero_sp", "Saccogl_ko", "Scutige_co", "Spizell_pu", "Sphaero_ar", "Strongy_pu",
                                                                   "Suberit_do", "Sycon_cili", "Sycon_raph", "Tethya_wil", "Trichop_ad", "Trichop_H4", "Tubifex_tu",
                                                                   "Xenotur_bo"),
                    "Nosenko2013.ribosomal_14615_smatrix.aa" = c("Acropor_hy", "Acropor_mi", "Acropor_pa", "Allomyc_ma", "Amoebid_pa", "Amphime_qu", "Anoplod_er",
                                                                 "Aplysia_ca", "Aphroca_va", "Asbesto_sp", "Astrosc_wi", "Batrach_de", "Beroe_sp", "Branchi_fl",
                                                                 "Capitel_sp", "Capsasp_ow", "Carteri_fo", "Ciona_inte", "Clytia_hem", "Corticium", "Crassos_gi",
                                                                 "Cratero_me", "Cyanea_cap", "Danio_reri", "Daphnia_pu", "Ephydat_mu", "Euprymn_sc", "Helobde_ro",
                                                                 "Heteroc_ca", "Hydra_magn", "Hydra_vulg", "Hydract_ec", "Ixodes_sca", "Leucett_ch", "Leucoso_sp",
                                                                 "Lubomir_ba", "Mertens_sp", "Metridi_se", "Mnemiop_le", "Molgula_te", "Monosig_br", "Monosig_ov",
                                                                 "Montast_fa", "Nasonia_vi", "Nematos_ve", "Oopsacas_m", "Oscarel_ca", "Oscarel_lo", "Oscarel_sp",
                                                                 "Pachydi_gl", "Pedicel_ce", "Pedicul_hu", "Petromy_ma", "Phycomy_bl", "Pleurob_pi", "Podocor_ca",
                                                                 "Porites", "Protero_sp", "Saccogl_ko", "Scutige_co", "Spizell_pu", "Sphaero_ar", "Strongy_pu",
                                                                 "Suberit_fu", "Suberit_do", "Sycon_raph", "Tethya_wil", "Trichop_ad", "Trichop_H4", "Tubifex_tu",
                                                                 "Xenotur_bo"))

# Create one list that collates datasets to contain all taxa for each clade 
all_taxa <- list("Bilateria" = sort(unique(unlist(lapply(all_datasets, function(d){d$Bilateria})))),
                 "Cnidaria" = sort(unique(unlist(lapply(all_datasets, function(d){d$Cnidaria})))),
                 "Placozoa" = sort(unique(unlist(lapply(all_datasets, function(d){d$Placozoa})))),
                 "Porifera" = sort(unique(unlist(lapply(all_datasets, function(d){d$Porifera})))),
                 "Ctenophora" = sort(unique(unlist(lapply(all_datasets, function(d){d$Ctenophora})))),
                 "Outgroup" = sort(unique(unlist(lapply(all_datasets, function(d){d$Outgroup})))),
                 "Sponges_Calcarea" = sort(unique(unlist(lapply(all_datasets, function(d){d$Sponges_Calcarea})))),
                 "Sponges_Homoscleromorpha" = sort(unique(unlist(lapply(all_datasets, function(d){d$Sponges_Homoscleromorpha})))),
                 "Sponges_Hexactinellida" = sort(unique(unlist(lapply(all_datasets, function(d){d$Sponges_Hexactinellida})))),
                 "Sponges_Demospongiae" = sort(unique(unlist(lapply(all_datasets, function(d){d$Sponges_Demospongiae})))) )



#### MODELS ####
models_list <- list("redmond2021" = c("C30", "C40", "C50", "C60", "EHO+G", "EX_EHO+G", "EX2+G", "EX3+G", 
                                      "JTT+C20", "JTT+C30", "JTT+C40", "JTT+C50", "JTT+C60", "JTT+G4", "JTTDCMut+G4", "JTTDCMut+I+G4", 
                                      "LG+C30", "LG+C40", "LG+C50", "LG+C60", "LG+CF4", "LG+F+G4", "LG+F+I+G4", "LG+G4", "LG+I+G4", "LG4M", 
                                      "mtZOA+F+G4", "mtZOA+G4", "mtZOA+I+G4", "PMB+G4", "rtREV+G4", "rtREV+I+G4", "UL2+G", "UL3+G", 
                                      "WAG+C20", "WAG+C30", "WAG+C40", "WAG+C50", "WAG+C60", "WAG+G4", "WAG+I+G4"),
                    "li2021" = c("C60+LG", "C60+Poisson", "C60+WAG", "CAT+F81", "CAT+GTR", "CAT+GTR+G", "CAT+Poisson",
                                 "CAT+Poisson+G", "CAT+WAG+F", "GTR+FO", "GTR+G+FO", "GTR20", "LG+G+F", "WAG",        
                                 "WAG+F", "WAG+G+F"),
                    "spreadsheet" = c("ModelFinder", "rtREV+G+F", "GTR+G4", 
                                      "WAG", "WAG+F", "WAG+CAT+G4", "WAG+C60", "WAG+C50", "WAG+C40", "WAG+C30", "WAG+C20", "WAG+C10",
                                      "LG", "LG+F", "LG+G", "LG+G4+F", "LG+C60", "LG+C50", "LG+C40", "LG+C30", "LG+C20", "LG+C10",
                                      "C20+LG+FO+R4", "C60+LG+FO+R4", "C60+LG+G+F", "LG+PMSF+G", 
                                      "Poisson+G", "Poisson+CAT", "Poisson+CAT+G", "JTT", "JTT+C60", "JTT+C50", "JTT+C40", "JTT+C30", "JTT+C20", "JTT+C10",
                                      "EX2+G", "EX3+G", "EX_EHO+G", "C60", "C50", "C40", "C30", "C20", "C10", 
                                      "UL3+G", "UL2+G", "mtZOA+G4", "LG4M"))
# Collate all models
all_models <- unlist(models_list)
# Sort and identify unique models
all_models <- sort(unique(all_models))


