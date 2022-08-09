## caitlinch/metazoan-mixtures/test_whelan2017_mixtures.R
# Caitlin Cherryh 2022

# Proof of concept for the mixture of trees method

### Step 1: Input parameters ###
# main_dir                <- path to caitlinch/metazoan-mixtures git repository
# gene_folder             <- path to folder containing fasta files for each gene in the Whelan 2017 dataset
# iqtree_path             <- path to IQ-Tree2 executable with mixtures of trees implementation
# constraint_tree_dir     <- folder to store constraint trees in
# number_parallel_threads <- number of cores to use for parallel processes

location = "local"
if (location == "local"){
  main_dir <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"
  gene_folder <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/01_Data_Whelan2017/genes/"
  short_gene_folder <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/01_Data_Whelan2017/genes_short/"
  iqtree_path <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/02_Software_IQ-Tree/IQ-Tree_2.2.0.3.tm.3/iqtree-2.2.0.3.tm.3-MacOSX/bin/iqtree2"
  constraint_tree_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_constraint_trees/"
  
  number_parallel_threads = 1
} else if (location == "soma"){
  main_dir <- "/data/caitlin/metazoan-mixtures/"
  gene_folder <- "/data/caitlin/metazoan-mixtures/data_whelan2017/genes/"
  iqtree_path <- "/data/caitlin/metazoan-mixtures/iqtree-2.2.0.3.tm.3-Linux/bin/iqtree2"
  constraint_tree_dir <- "/data/caitlin/metazoan-mixtures/constraint_trees/"
  
  number_parallel_threads = 20
}

assemble_constraint_trees <- FALSE
estimate_constraint_trees <- TRUE



### Step 2: Prepare analysis ###
# Source function files
source(paste0(main_dir, "code/func_constraint_trees.R"))

# Open packages
library(parallel)
library(ape)
library(phangorn)

# Create folders if necessary
if (dir.exists(constraint_tree_dir) == FALSE){dir.create(constraint_tree_dir)}



### Step 3: Prepare constraint trees ###
if (assemble_constraint_trees == TRUE){
  ## For Whelan2017 data:
  # Set dataset name
  dataset = "Whelan2017"
  
  # Create folder for each dataset inside the constraint tree folder
  dataset_constraint_tree_dir <- paste0(constraint_tree_dir, dataset, "/")
  if (dir.exists(dataset_constraint_tree_dir) == FALSE){dir.create(dataset_constraint_tree_dir)}
  
  # All taxa
  # Identify which taxa are in which clades
  bilateria_taxa = c("Homo_sapiens", "Strongylocentrotus_purpatus", "Hemithris_psittacea", "Capitella_teleta", "Drosophila_melanogaster",
                     "Daphnia_pulex")
  cnidaria_taxa = c("Hydra_vulgaris", "Bolocera_tuediae", "Aiptasia_pallida", "Hormathia_digitata", "Nematostella_vectensis", "Acropora_digitifera",
                    "Eunicella_verrucosa", "Hydra_viridissima", "Hydra_oligactis", "Physalia_physalia", "Abylopsis_tetragona",
                    "Craseo_lathetica", "Nanomia_bijuga", "Agalma_elegans", "Periphyla_periphyla")
  placozoa_taxa = c("Trichoplax_adhaerens")
  porifera_taxa = c("Cliona_varians", "Sycon_coactum", "Sycon_ciliatum", "Corticium_candelabrum", "Oscarella_carmela", 
                    "Hyalonema_populiferum", "Aphrocallistes_vastus", "Rossella_fibulata", "Sympagella_nux", "Ircinia_fasciculata",
                    "Chondrilla_nucula", "Amphimedon_queenslandica", "Petrosia_ficiformis", "Spongilla_lacustris", 
                    "Pseudospongosorites_suberitoides", "Mycale_phylophylla", "Latrunculia_apicalis", "Crella_elegans",
                    "Kirkpatrickia_variolosa")
  ctenophora_taxa = c("Euplokamis_dunlapae", "Vallicula_sp", "Coeloplana_astericola", "Hormiphora_californica", "Hormiphora_palmata",
                      "Pleurobrachia_pileus", "Pleurobrachia_bachei", "Pleurobrachia_sp_South_Carolina_USA", "Cydippida_sp_Maryland_USA",
                      "Callianira_Antarctica", "Mertensiidae_sp_Antarctica", "Mertensiidae_sp_Washington_USA", "Cydippida_sp",
                      "Dryodora_glandiformis", "Lobatolampea_tetragona", "Beroe_abyssicola", "Beroe_sp_Antarctica", "Beroe_ovata",
                      "Beroe_sp_Queensland_Australia", "Beroe_forskalii", "Ocyropsis_sp_Bimini_Bahamas", "Ocyropsis_crystallina",
                      "Ocyropsis_sp_Florida_USA", "Bolinopsis_infundibulum", "Mnemiopsis_leidyi", "Bolinopsis_ashleyi", 
                      "Lobata_sp_Punta_Arenas_Argentina", "Eurhamphaea_vexilligera", "Cestum_veneris", "Ctenophora_sp_Florida_USA")
  outgroup_taxa = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis")
  # Break sponges into two groups to allow for paraphyly - based on Whelan et al (2017) main figure
  sponges_1_taxa = c("Sycon_coactum", "Sycon_ciliatum", "Oscarella_carmela", "Corticium_candelabrum")
  sponges_2.1_taxa = c("Hyalonema_populiferum", "Sympagella_nux", "Rossella_fibulata", "Aphrocallistes_vastus")
  sponges_2.2_taxa = c("Ircinia_fasciculata", "Chondrilla_nucula", "Spongilla_lacustris", "Cliona_varians",
                       "Pseudospongosorites_suberitoides", "Mycale_phylophylla", "Latrunculia_apicalis", 
                       "Kirkpatrickia_variolosa", "Crella_elegans", "Petrosia_ficiformis", "Amphimedon_queenslandica")
  sponges_2_taxa = c(sponges_2.1_taxa, sponges_2.2_taxa)
  
  # Break sponges into groups according to classification
  # Note: Homoscleromorpha was thought to belong to the Demospongiae, but it is phylogenetically well-separated.
  #       See: Bergquist PR (1978). Sponges. London: Hutchinson. ISBN 978-0-520-03658-1.
  ## Porifera topologies:
  # Borchiellini et al (2001) found Porifera paraphyletic, the Calcarea being more related to monophyletic Eumetazoa than to the 
  #     siliceous sponges (Demospongiae, Hexactinellida)
  #     QUOTE: "Hexactinellida appears to be the sister-group of a Demospongiae–Calcarea–Cnidaria–Ctenophora–Placozoa clade both in 
  #      neighbour-joining and in maximum parsimony"
  # Medina et al (2001) found Porifera paraphyletic. Either Calcarea as sister to all animals, then (Hexactinellida, Demospongiae) or
  #     vice versa (i.e. (Hexactinellida, Demospongiae) first then Calcarea) 
  #     Note: placement of Ctenophora different in their two trees (one from SSU rRNA, one from LSU rRNA)
  # Sperling et al (2007) find paraphyletic sponges. Demosponges monophyletic and sister to all other animals. Then branches off 
  #     Calcarea, then Homoscleromorpha.
  # Dohrmann et al (2008) finds monophyletic Porifera, but the clade consists of two sister clades 
  #     (Calcarea, Homoscleromorpha), (Hexactinellida, Demospongiae))
  # Sperling et al (2009) has paraphlyetic Porifera. One clade combines (Hexactinellida, Demospongiae) as sister to all animals. 
  #     Then branches off Calcarea, then Homoscleromorpha.
  # Whelan et al (2017) tree has monophyletic Porifera, but the clade consists of two sister clades 
  #     (Calcarea, Homoscleromorpha), (Hexactinellida, Demospongiae))
  # Sperling et al (2010) find paraphlyetic Porifera. Demosponges are monophyletic, and that hexactinellids are their sister group
  #     (together forming the Silicea as sister to all animals). Then branches offCalcarea, then Homoscleromorpha.
  sponges_calcarea_taxa = c("Sycon_coactum", "Sycon_ciliatum")
  sponges_homoscleromorpha_taxa = c("Oscarella_carmela", "Corticium_candelabrum")
  sponges_hexactinellida_taxa = c("Hyalonema_populiferum", "Sympagella_nux", "Rossella_fibulata", "Aphrocallistes_vastus")
  sponges_demospongiae_taxa = c("Ircinia_fasciculata", "Chondrilla_nucula", "Spongilla_lacustris", "Cliona_varians", 
                                "Pseudospongosorites_suberitoides", "Mycale_phylophylla", "Latrunculia_apicalis",
                                "Kirkpatrickia_variolosa", "Crella_elegans", "Petrosia_ficiformis", "Amphimedon_queenslandica")
  # Hypotheses:
  #   1. Ctenophora-sister
  #   2. Porifera-sister
  #   3. Porifera+Ctenophora-sister
  #   4. Paraphyletic sponges, Porifera-sister
  #   5. Paraphyletic sponges, Ctenophora-sister
  # Uninvestigated hypotheses:
  #   1. Placozoa-sister
  
  ## Hypothesis 1: Ctenophora-sister
  # Tree: (outgroup_taxa, (ctenophora_taxa, (porifera_taxa, (placozoa_taxa, cnidaria_taxa, bilateria_taxa))))
  # Construct constraint tree
  constraint_tree_1 <- paste0("((", 
                              paste(outgroup_taxa, collapse = ", "), 
                              "),((", 
                              paste(ctenophora_taxa, collapse = ", "), 
                              "),((", 
                              paste(c(porifera_taxa), collapse = ", "),
                              "), (",
                              paste(c(placozoa_taxa, cnidaria_taxa, bilateria_taxa), collapse = ", "),
                              "))));")
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, dataset, "_constraint_tree_", "1", ".nex")
  write(constraint_tree_1, file = constraint_tree_file_name)
  
  
  ## Hypothesis 2: Porifera-sister
  # Tree: (outgroup_taxa, (porifera_taxa, (ctenophora_taxa, (placozoa_taxa, cnidaria_taxa, bilateria_taxa))))
  # Construct constraint tree
  constraint_tree_2 <- paste0("((", 
                              paste(outgroup_taxa, collapse = ", "), 
                              "),((", 
                              paste(porifera_taxa, collapse = ", "), 
                              "),((", 
                              paste(c(ctenophora_taxa), collapse = ", "),
                              "),(",
                              paste(c(placozoa_taxa, cnidaria_taxa, bilateria_taxa), collapse = ", "), 
                              "))));")
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, dataset, "_constraint_tree_", "2", ".nex")
  write(constraint_tree_2, file = constraint_tree_file_name)
  
  
  ## Hypothesis 3: Porifera+Ctenophora-sister
  # Tree: (outgroup_taxa, ((porifera_taxa, ctenophora_taxa), (placozoa_taxa, (cnidaria_taxa, bilateria_taxa))))
  # Construct constraint tree
  constraint_tree_3 <- paste0("((", 
                              paste(outgroup_taxa, collapse = ", "), 
                              "),((", 
                              paste(c(porifera_taxa, ctenophora_taxa), collapse = ", "), 
                              "),(", 
                              paste(c(placozoa_taxa), collapse = ", "),
                              ", (",
                              paste(c(cnidaria_taxa, bilateria_taxa), collapse = ", "), 
                              "))));")
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, dataset, "_constraint_tree_", "3", ".nex")
  write(constraint_tree_3, file = constraint_tree_file_name)
  
  ## Hypothesis 4: Paraphyletic sponges, Porifera-sister
  # Tree: (outgroup_taxa, (sponges_1_taxa, (sponges_2_taxa, (ctenophora_taxa, (placozoa_taxa, cnidaria_taxa, bilateria_taxa)))))
  # Construct constraint tree
  constraint_tree_4 <- paste0("((", 
                              paste(outgroup_taxa, collapse = ", ") ,
                              ") ,((", 
                              paste(sponges_1_taxa, collapse = ", "), 
                              "), ((", 
                              paste(sponges_2_taxa, collapse = ", "), 
                              "), ((", 
                              paste(c(ctenophora_taxa), collapse = ", "),
                              "), (",
                              paste(c(placozoa_taxa, cnidaria_taxa, bilateria_taxa), collapse = ", "),
                              ")))));")
  write(constraint_tree_4, file = constraint_tree_file_name)
  
  ## Hypothesis 5: Paraphyletic sponges, Ctenophora-sister
  # Tree: (outgroup_taxa, (ctenophora_taxa, (sponges_1_taxa, (sponges_2_taxa, (placozoa_taxa, cnidaria_taxa, bilateria_taxa)))))
  # Construct constraint tree
  constraint_tree_5 <- paste0("((", 
                              paste(outgroup_taxa, collapse = ", "),
                              ") ,((",
                              paste(ctenophora_taxa, collapse = ", "),
                              "), ((", 
                              paste(sponges_1_taxa, collapse = ", "),
                              "), ((", 
                              paste(c(sponges_2_taxa), collapse = ", "),
                              "), (",
                              paste(c(placozoa_taxa, cnidaria_taxa, bilateria_taxa), collapse = ", "),
                              ")))));")
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, dataset, "_constraint_tree_", "5", ".nex")
  write(constraint_tree_5, file = constraint_tree_file_name)
  
  # Assemble dataframe of information about the constraint trees
  constraint_df <- data.frame(constraint_tree_id = 1:5,
                              constraint_tree_paths = paste0(dataset_constraint_tree_dir, dataset, "_constraint_tree_", 1:5, ".nex"),
                              constraint_prefixes = paste0(dataset, "_ConstraintTree", 1:5),
                              alignment_path = gene_folder,
                              model = "Q.insect+R8",
                              iqtree_path = iqtree_path,
                              constraint_trees = c(constraint_tree_1, constraint_tree_2, constraint_tree_3, 
                                                   constraint_tree_4, constraint_tree_5),
                              num_threads = number_parallel_threads)
  
  # Write dataset of information about constraint trees
  constraint_df_path <- paste0(dataset_constraint_tree_dir, dataset, "_constraint_tree_parameters.csv")
  write.csv(constraint_df, constraint_df_path, row.names = FALSE)
  
  ## Small trees test
  # Identify which taxa are in which clades
  bilateria_taxa = c("Homo_sapiens", "Drosophila_melanogaster")
  cnidaria_taxa = c("Hydra_vulgaris", "Hydra_viridissima")
  placozoa_taxa = c("Trichoplax_adhaerens")
  porifera_taxa = c("Sycon_coactum", "Sycon_ciliatum", "Sympagella_nux", "Pseudospongosorites_suberitoides")
  ctenophora_taxa = c("Mnemiopsis_leidyi", "Ctenophora_sp_Florida_USA")
  outgroup_taxa = c("Monosiga_ovata", "Monosiga_brevicolis")
  # Break sponges into two groups to allow for paraphyly - based on Whelan et al (2017) main figure
  sponges_1_taxa = c("Sycon_coactum", "Sycon_ciliatum")
  sponges_2.1_taxa = c("Sympagella_nux")
  sponges_2.2_taxa = c("Pseudospongosorites_suberitoides")
  sponges_2_taxa = c(sponges_2.1_taxa, sponges_2.2_taxa)
  
  # Hypotheses:
  #   1. Ctenophora-sister
  #   2. Porifera-sister
  #   3. Porifera+Ctenophora-sister
  #   4. Paraphyletic sponges, Porifera-sister
  #   5. Paraphyletic sponges, Ctenophora-sister
  # Uninvestigated hypotheses:
  #   1. Placozoa-sister
  
  ## Hypothesis 1: Ctenophora-sister
  # Tree: (outgroup_taxa, (ctenophora_taxa, (porifera_taxa, (placozoa_taxa, cnidaria_taxa, bilateria_taxa))))
  # Construct constraint tree
  constraint_tree_1 <- paste0("((", 
                              paste(outgroup_taxa, collapse = ", "), 
                              "),((", 
                              paste(ctenophora_taxa, collapse = ", "), 
                              "),((", 
                              paste(c(porifera_taxa), collapse = ", "),
                              "), (",
                              paste(c(placozoa_taxa, cnidaria_taxa, bilateria_taxa), collapse = ", "),
                              "))));")
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, dataset, "_small_constraint_tree_", "1", ".nex")
  write(constraint_tree_1, file = constraint_tree_file_name)
  
  
  ## Hypothesis 2: Porifera-sister
  # Tree: (outgroup_taxa, (porifera_taxa, (ctenophora_taxa, (placozoa_taxa, cnidaria_taxa, bilateria_taxa))))
  # Construct constraint tree
  constraint_tree_2 <- paste0("((", 
                              paste(outgroup_taxa, collapse = ", "), 
                              "),((", 
                              paste(porifera_taxa, collapse = ", "), 
                              "),((", 
                              paste(c(ctenophora_taxa), collapse = ", "),
                              "),(",
                              paste(c(placozoa_taxa, cnidaria_taxa, bilateria_taxa), collapse = ", "), 
                              "))));")
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, dataset, "_small_constraint_tree_", "2", ".nex")
  write(constraint_tree_2, file = constraint_tree_file_name)
  
  
  ## Hypothesis 3: Porifera+Ctenophora-sister
  # Tree: (outgroup_taxa, ((porifera_taxa, ctenophora_taxa), (placozoa_taxa, (cnidaria_taxa, bilateria_taxa))))
  # Construct constraint tree
  constraint_tree_3 <- paste0("((", 
                              paste(outgroup_taxa, collapse = ", "), 
                              "),((", 
                              paste(c(porifera_taxa, ctenophora_taxa), collapse = ", "), 
                              "),(", 
                              paste(c(placozoa_taxa), collapse = ", "),
                              ", (",
                              paste(c(cnidaria_taxa, bilateria_taxa), collapse = ", "), 
                              "))));")
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, dataset, "_small_constraint_tree_", "3", ".nex")
  write(constraint_tree_3, file = constraint_tree_file_name)
  
  ## Hypothesis 4: Paraphyletic sponges, Porifera-sister
  # Tree: (outgroup_taxa, (sponges_1_taxa, (sponges_2_taxa, (ctenophora_taxa, (placozoa_taxa, cnidaria_taxa, bilateria_taxa)))))
  # Construct constraint tree
  constraint_tree_4 <- paste0("((", 
                              paste(outgroup_taxa, collapse = ", ") ,
                              ") ,((", 
                              paste(sponges_1_taxa, collapse = ", "), 
                              "), ((", 
                              paste(sponges_2_taxa, collapse = ", "), 
                              "), ((", 
                              paste(c(ctenophora_taxa), collapse = ", "),
                              "), (",
                              
                              paste(c(placozoa_taxa, cnidaria_taxa, bilateria_taxa), collapse = ", "),
                              ")))));")
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, dataset, "_small_constraint_tree_", "4", ".nex")
  write(constraint_tree_4, file = constraint_tree_file_name)
  
  ## Hypothesis 5: Paraphyletic sponges, Ctenophora-sister
  # Tree: (outgroup_taxa, (ctenophora_taxa, (sponges_1_taxa, (sponges_2_taxa, (placozoa_taxa, cnidaria_taxa, bilateria_taxa)))))
  # Construct constraint tree
  constraint_tree_5 <- paste0("((", 
                              paste(outgroup_taxa, collapse = ", "),
                              ") ,((",
                              paste(ctenophora_taxa, collapse = ", "),
                              "), ((", 
                              paste(sponges_1_taxa, collapse = ", "),
                              "), ((", 
                              paste(c(sponges_2_taxa), collapse = ", "),
                              "), (",
                              paste(c(placozoa_taxa, cnidaria_taxa, bilateria_taxa), collapse = ", "),
                              ")))));")
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, dataset, "_small_constraint_tree_", "5", ".nex")
  write(constraint_tree_5, file = constraint_tree_file_name)
  
  # Assemble dataframe of information about the constraint trees
  constraint_df <- data.frame(constraint_tree_id = 1:5,
                              constraint_tree_paths = paste0(dataset_constraint_tree_dir, dataset, "_small_constraint_tree_", 1:5, ".nex"),
                              constraint_prefixes = paste0(dataset, "_SmallConstraintTree", 1:5),
                              alignment_path = short_gene_folder,
                              model = c("Q.insect+R8", "Q.insect+R8", "Q.insect+R8", "Q.insect+R8", "Q.insect+R8"),
                              iqtree_path = iqtree_path,
                              constraint_trees = c(constraint_tree_1, constraint_tree_2, constraint_tree_3, 
                                                   constraint_tree_4, constraint_tree_5),
                              num_threads = "AUTO")
  
  # Write dataset of information about constraint trees
  constraint_df_path <- paste0(dataset_constraint_tree_dir, dataset, "_small_constraint_tree_parameters.csv")
  write.csv(constraint_df, constraint_df_path, row.names = FALSE)
  
}



### Step 4: Estimate trees with constraint trees ###
if (estimate_constraint_trees == TRUE){
  ## For Whelan2017 data:
  # Set dataset name
  dataset = "Whelan2017"
  
  # Create folder for each dataset inside the constraint tree folder
  dataset_constraint_tree_dir <- paste0(constraint_tree_dir, dataset, "/")
  if (dir.exists(dataset_constraint_tree_dir) == FALSE){dir.create(dataset_constraint_tree_dir)}
  
  
  # Set working directory to dataset_constraint_tree_dir so IQ-Tree output is saved with the constraint trees
  setwd(dataset_constraint_tree_dir)
  
  # For trees with all 76 taxa
  estimate_trees_df <- read.csv(paste0(dataset_constraint_tree_dir, dataset, "_constraint_tree_parameters.csv"))
  # Estimate an ML tree in IQ-Tree for each constraint tree
  lapply(1:nrow(estimate_trees_df), apply.one.constraint.tree, estimate_trees_df)
}



### Step 5: Compare trees
# Open Whelan 2017 constraint trees
all_dir_files <- list.files(constraint_tree_dir, recursive = TRUE)
whelan2017_files <- grep("Whelan2017", all_dir_files, value = TRUE)
constraint_tree_files <- grep(".treefile", whelan2017_files, value = TRUE)
small_constraint_trees <- grep("Small", constraint_tree_files, value = TRUE)
ctrees <- read.tree(text = unlist(lapply(paste0(constraint_tree_dir, small_constraint_trees), readLines))) #open small constraint trees
rooted_ctrees <- root(ctrees, c("Monosiga_ovata", "Monosiga_brevicolis"))





