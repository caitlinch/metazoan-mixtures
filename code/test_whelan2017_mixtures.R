## caitlinch/metazoan-mixtures/test_whelan2017_mixtures.R
# Caitlin Cherryh 2022

# Proof of concept for the mixture of trees method

### Step 1: Input parameters ###
# gene_folder <- path to folder containing fasta files for each gene in the Whelan 2017 dataset
# iqtree_path <- path to IQ-Tree2 executable with mixtures of trees implementation

gene_folder <- "/data/caitlin/metazoan-mixtures/data_whelan2017/genes/"
iqtree_path <- "/data/caitlin/metazoan-mixtures/iqtree-2.2.0.3.tm.3-Linux/bin/iqtree2"



### Step 2: Prepare analysis###
# Open packages
library(ape)

# Identify which taxa are in which clades
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



### Step 3: Prepare tree topology hypotheses ###
# Hypotheses:
#   1. Ctenophora-sister
#   2. Porifera-sister
#   3. Porifera+Ctenophora-sister

