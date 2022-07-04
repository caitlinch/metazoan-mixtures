# Give location for genes
gene_folder <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/01_Data_Whelan2017/genes"
# Give location of new IQ-Tree version
iqtree_path <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/02_Software_IQ-Tree/IQ-Tree_2.2.0.3.tm.3/iqtree-2.2.0.3.tm.3-MacOSX/bin/iqtree2"
# Give location of trees
tree_file <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_trees/Whelan2017/Species_tree/partitions.nex.treefile"
# Output location for clades
clade_tree_folder <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_trees/Whelan2017/Clades/"

# Open packages
library(ape)

# Open tree
tree <- read.tree(tree_file)
# List taxa in each clade
bilateria_taxa = c("Homo_sapiens", "Strongylocentrotus_purpatus", "Hemithris_psittacea", "Capitella_teleta", "Drosophila_melanogaster",
              "Daphnia_pulex")
cnidaria_taxa = c("Hydra_vulgaris", "Bolocera_tuediae", "Aiptasia_pallida", "Hormathia_digitata", "Nematostella_vectensis", "Acropora_digitifera",
             "Eunicella_verrucosa", "Hydra_viridissima", "Hydra_oligactis", "Physalia_physalia", "Abylopsis_tetragona",
             "Craseo_lathetica", "Nanomia_bijuga", "Agalma_elegans", "Periphyla_periphyla")
placozoa_taxa = c("Trichoplax_adhaerens")
porifera_taxa = c("Cliona_varians", "Sycon_coactum", "Sycon_ciliatum", "Corticium_candelabrum", "Oscarella_carmela", "Hyalonema_populiferum",
             "Aphrocallistes_vastus", "Rossella_fibulata", "Sympagella_nux", "Ircinia_fasciculata", "Chondrilla_nucula", "Amphimedon_queenslandica",
             "Petrosia_ficiformis", "Spongilla_lacustris", "Pseudospongosorites_suberitoides", "Mycale_phylophylla", "Latrunculia_apicalis", "Crella_elegans",
             "Kirkpatrickia_variolosa")
ctenophora_taxa = c("Euplokamis_dunlapae", "Vallicula_sp", "Coeloplana_astericola", "Hormiphora_californica", "Hormiphora_palmata",
               "Pleurobrachia_pileus", "Pleurobrachia_bachei", "Pleurobrachia_sp_South_Carolina_USA", "Cydippida_sp_Maryland_USA",
               "Callianira_Antarctica", "Mertensiidae_sp_Antarctica", "Mertensiidae_sp_Washington_USA", "Cydippida_sp",
               "Dryodora_glandiformis", "Lobatolampea_tetragona", "Beroe_abyssicola", "Beroe_sp_Antarctica", "Beroe_ovata",
               "Beroe_sp_Queensland_Australia", "Beroe_forskalii", "Ocyropsis_sp_Bimini_Bahamas", "Ocyropsis_crystallina",
               "Ocyropsis_sp_Florida_USA", "Bolinopsis_infundibulum", "Mnemiopsis_leidyi", "Bolinopsis_ashleyi", 
               "Lobata_sp_Punta_Arenas_Argentina", "Eurhamphaea_vexilligera", "Cestum_veneris", "Ctenophora_sp_Florida_USA")
outgroup_taxa = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis")
# Drop tips to make a taxa of each clade
bilateria_clade <- keep.tip(tree, bilateria_taxa)
cnidaria_clade <- keep.tip(tree, cnidaria_taxa)
placozoa_clade <- keep.tip(tree, placozoa_taxa)
porifera_clade <- keep.tip(tree, porifera_taxa)
ctenophora_clade <- keep.tip(tree, ctenophora_taxa)
outgroup_clade <- keep.tip(tree, outgroup_taxa)
# Remove node values
bilateria_clade$node.label <- NULL
cnidaria_clade$node.label <- NULL
placozoa_clade$node.label <- NULL
porifera_clade$node.label <- NULL
ctenophora_clade$node.label <- NULL
outgroup_clade$node.label <- NULL
# Remove node values
bilateria_clade$edge.length <- NULL
cnidaria_clade$edge.length <- NULL
placozoa_clade$edge.length <- NULL
porifera_clade$edge.length <- NULL
ctenophora_clade$edge.length <- NULL
outgroup_clade$edge.length <- NULL
# Save topology only of each clade
write.tree(bilateria_clade, file = paste0(clade_tree_folder, "bilateria_clade.tre"))
write.tree(cnidaria_clade, file = paste0(clade_tree_folder, "cnidaria_clade.tre"))
write.tree(placozoa_clade, file = paste0(clade_tree_folder, "placozoa_clade.tre"))
write.tree(porifera_clade, file = paste0(clade_tree_folder, "porifera_clade.tre"))
write.tree(ctenophora_clade, file = paste0(clade_tree_folder, "ctenophora_clade.tre"))
write.tree(outgroup_clade, file = paste0(clade_tree_folder, "outgroup_clade.tre"))
