## caitlinch/metazoan-mixtures/code/03_data_analysis_and_plotting.R
# This script performs data analysis and creates plots of methods and results
# Caitlin Cherryh 2023

## This script:
# 1. Plots figures and trees for the manuscript


#### 1. Input parameters ####
## Specify parameters:
# plot_dir          <- Directory for IQ-Tree output (trees and tree mixtures)
# results_dir         <- Directory for results and plots
# repo_dir            <- Location of caitlinch/metazoan-mixtures github repository

plot_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/05_plotting/"
results_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/01_output_files/"
tree_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/02_maximum_likelihood_trees/02_renamed_trees/"
repo_dir <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"



#### 2. Prepare functions, variables and packages ####
# Open packages
library(ape)
library(ggplot2)
library(ggtree)
library(patchwork)
library(reshape2)

# Open function files
source(paste0(repo_dir,"code/func_plotting.R"))
source(paste0(repo_dir,"code/func_data_processing.R"))



#### 3. Plot figures for introduction and methods sections of manuscript ####
### Plotting the alternative phylogenetic hypotheses ###
# Open the trees
trees_file <- paste0(repo_dir, "trees/alternative_phylogenetic_hypotheses.nex")
trees <- read.tree(trees_file)
image_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/99_RSB_HDR_conference_2023/pictures/"

# Create a dataframe for labelling
metazoan_clade_labs <- data.frame("clade" = c("Bilateria", "Cnidaria", "Ctenophora", "Porifera", "Calcarea", "Demospongiae", "Hexactinellida", "Homoscleromorpha", "Outgroup"),
                                  "color" = c("A", "B", "C", rep("D", 5), "E"))
metazoan_clade_labs$lab <- metazoan_clade_labs$clade
metazoan_clade_labs$html_lab <- paste0("<b style='color:", metazoan_clade_labs$color, "'>", metazoan_clade_labs$clade, "</b>")

### Black and White plots ###
# Plot one tree at a time
hyp1_plot <- bw.monophyletic.clades.plot(trees[[1]])
hyp2_plot <- bw.monophyletic.clades.plot(trees[[2]])
hyp3_plot <- bw.monophyletic.clades.plot(trees[[3]])
hyp4_plot <- bw.paraphyletic.clades.plot(trees[[4]], label_nodes = c(12, 14))
hyp5_plot <- bw.paraphyletic.clades.plot(trees[[5]], label_nodes = c(11, 13))
# Combine the four hypothesis trees into one 
patchwork_hyps <- (hyp1_plot | hyp2_plot | hyp3_plot)/(hyp4_plot | hyp5_plot) + 
  plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 20))
# Export hypothesis plots 
hypothesis_plot_file <- paste0(plot_dir, "hypothesis_tree_example_plot.png")
png(filename = hypothesis_plot_file, width = 1090, height = 723, units = "px", pointsize = 12, bg = "white")
patchwork_hyps
dev.off()


### Colour plots ###
# Colourblind friendly palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
metazoan_palette <- c(A = "#CC79A7", B = "#009E73", C = "#56B4E9", D = "#E69F00", E = "#999999")
tree2_palette <- c("#F0F921FF", "#0D0887FF")
tree5_palette <- c("#e66101", "#fdb863", "#f7f7f7", "#b2abd2", "#5e3c99")
tree5_cividis <- c("#FDE725FF", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF")
tree2_cividis <- c(tree5_cividis[1], tree5_cividis[5])
tree2_tonal <- c("#bdd7e7", "#2171b5")

# Plot one tree at a time
hyp1_plot <- color.clades.plot(trees[[1]], tip_labels = metazoan_clade_labs, color_palette = metazoan_palette, 
                               save.plot = TRUE, output_directory = plot_dir, output_id = "hypothesis_tree_1_plot_color")
hyp2_plot <- color.clades.plot(trees[[2]], tip_labels = metazoan_clade_labs, color_palette = metazoan_palette, 
                               save.plot = TRUE, output_directory = plot_dir, output_id = "hypothesis_tree_2_plot_color")
hyp3_plot <- color.clades.plot(trees[[3]], tip_labels = metazoan_clade_labs, color_palette = metazoan_palette, 
                               save.plot = TRUE, output_directory = plot_dir, output_id = "hypothesis_tree_3_plot_color")
hyp4_plot <- color.clades.plot(trees[[4]], tip_labels = metazoan_clade_labs, color_palette = metazoan_palette, 
                               save.plot = TRUE, output_directory = plot_dir, output_id = "hypothesis_tree_4_plot_color")
hyp5_plot <- color.clades.plot(trees[[5]], tip_labels = metazoan_clade_labs, color_palette = metazoan_palette, 
                               save.plot = TRUE, output_directory = plot_dir, output_id = "hypothesis_tree_5_plot_color")



#### 4. Exploratory plots ####
# List all files in the results directory
all_files <- list.files(results_dir)
# Open the summary alignment details and the summary topology results
al_df <- read.csv(paste0(results_dir, grep("summary_alignment_details", all_files, value = T)), stringsAsFactors = FALSE)
al_df$matrix_name[which(al_df$dataset == "Whelan2015")] <- "Dataset10"
topo_df <-read.csv(paste0(results_dir, grep("summary_ML_tree_topology", all_files, value = T)), stringsAsFactors = FALSE)
topo_df$ID <- paste0(topo_df$dataset, ".", topo_df$matrix_name)
al_df$ID <- paste0(al_df$dataset, ".", al_df$matrix_name)
al_df_row_order <- match(topo_df$ID, al_df$ID)[which(!is.na(match(topo_df$ID, al_df$ID)))]
al_df <- al_df[al_df_row_order, ]
# Make a df for plotting
al_df$percent_CTEN_sister <- topo_df$percent_CTEN_sister
al_df$percent_PORI_sister <- topo_df$percent_PORI_sister
al_df$percent_CTEN_PORI_sister <- topo_df$`percent_CTEN.PORI_sister`
al_df$percent_Radiata_sister <- topo_df$percent_Radiata_sister
al_df$percent_PORI_monophyletic <- topo_df$percent_PORI_monophyletic
al_df$percent_PORI_paraphyletic <- topo_df$percent_PORI_paraphyletic
al_df$percent_PORI_one_taxon <- topo_df$percent_PORI_one_taxon
al_df$percent_CTEN_CNID_monophyletic <- topo_df$percent_CTEN.CNID_monophyletic
al_df$percent_CTEN_CNID_not_monophyletic <- topo_df$percent_CTEN.CNID_not_monophyletic
al_df$ID <- c("Dunn2008", "Philippe2009", "Pick2010", "Philippe2011", "Nosenko2013 non-ribo",
              "Nosenko2013 ribo", "Ryan2013", "Moroz2014", "Borowiec2015", "Chang2015",
              "Whelan2015", "Whelan2017", "Laumer2018", "Laumer2019")
al_df$best_model <- c("PMSF_C60", "PMSF_C60", "PMSF_C60", "PMSF_C60", "PMSF_C60",
                      "PMSF_C60", "PMSF_C60", "PMSF_C60", "PMSF_C60", "PMSF_LG_C60",
                      "PMSF_C60", "PMSF_C60", "PMSF_C60", "PMSF_C60")
# Extract branch a and branch b lengths
al_df$ctenophora_branch_length <- unlist(lapply(1:nrow(al_df), extract.branch.length.wrapper, alignment_df = al_df, tree_directory = tree_dir, clade = "Ctenophora"))
al_df$porifera_branch_length <- unlist(lapply(1:nrow(al_df), extract.branch.length.wrapper, alignment_df = al_df, tree_directory = tree_dir, clade = "Porifera"))



### Plot number of sites/number of informative sites against proportion of trees with each topology ###
### Plot 1: Percent of trees with Ctenophora sister against number of sites ###
plot_df <- melt(al_df, 
                id.vars = c("ID", "best_model", "num_taxa", "num_sites", "number_constant_sites", "proportion_constant_sites", 
                            "number_invariant_sites", "proportion_invariant_sites", "number_informative_sites",
                            "proportion_informative_sites"),
                measure.vars = c("percent_CTEN_sister") )
# Number of sites
ggplot(data = plot_df, aes(x = num_sites, y = value)) +
  geom_point() +
  geom_smooth() +
  scale_y_continuous(name ="Percentage of trees with Ctenophora-sister", breaks = seq(0,100,10), limits = c(0,110)) +
  scale_x_continuous(name = "Number of sites", breaks = seq(0,90000,15000), minor_breaks = seq(0,90000,5000)) +
  theme_bw()
# Number of constant sites
ggplot(data = plot_df, aes(x = proportion_constant_sites, y = value)) +
  geom_point() +
  scale_y_continuous(name ="Percentage of trees with Ctenophora-sister", breaks = seq(0,100,10), limits = c(0,110)) +
  scale_x_continuous(name = "Proportion of constant sites", breaks = seq(0.1,0.35,0.05), limits = c(0.1,0.35)) +
  theme_bw()
# Number of informative sites
ggplot(data = plot_df, aes(x = proportion_informative_sites, y = value)) +
  geom_point() +
  scale_y_continuous(name ="Percentage of trees with Ctenophora-sister", breaks = seq(0,100,10), limits = c(0,110)) +
  scale_x_continuous(name = "Proportion of informative sites", breaks = seq(0.5,0.8,0.05)) +
  theme_bw()

### Plot 2: Percent of trees with monophyletic Porifer against number of sites ###
plot_df <- melt(al_df, 
                id.vars = c("ID", "num_taxa", "num_sites", "number_constant_sites", "proportion_constant_sites", 
                            "number_invariant_sites", "proportion_invariant_sites", "number_informative_sites",
                            "proportion_informative_sites"),
                measure.vars = c("percent_PORI_monophyletic") )
# Number of sites
ggplot(data = plot_df, aes(x = num_sites, y = value)) +
  geom_point() +
  scale_y_continuous(name ="Percentage of trees with Ctenophora-sister", breaks = seq(0,100,10), limits = c(0,110)) +
  scale_x_continuous(name = "Number of sites", breaks = seq(0,90000,15000), minor_breaks = seq(0,90000,5000)) +
  theme_bw()
# Number of constant sites
ggplot(data = plot_df, aes(x = proportion_constant_sites, y = value)) +
  geom_point() +
  scale_y_continuous(name ="Percentage of trees with Ctenophora-sister", breaks = seq(0,100,10), limits = c(0,110)) +
  scale_x_continuous(name = "Proportion of constant sites", breaks = seq(0.1,0.35,0.05), limits = c(0.1,0.35)) +
  theme_bw()
# Number of informative sites
ggplot(data = plot_df, aes(x = proportion_informative_sites, y = value)) +
  geom_point() +
  scale_y_continuous(name ="Percentage of trees with Ctenophora-sister", breaks = seq(0,100,10), limits = c(0,110)) +
  scale_x_continuous(name = "Proportion of informative sites", breaks = seq(0.5,0.8,0.05)) +
  theme_bw()

### Plot 3: Percent of trees with Cten/Cnid paraphyletic against number of sites ###
plot_df <- melt(al_df, 
                id.vars = c("ID", "num_taxa", "num_sites", "number_constant_sites", "proportion_constant_sites", 
                            "number_invariant_sites", "proportion_invariant_sites", "number_informative_sites",
                            "proportion_informative_sites"),
                measure.vars = c("percent_CTEN_CNID_not_monophyletic") )
# Number of sites
ggplot(data = plot_df, aes(x = num_sites, y = value)) +
  geom_point() +
  scale_y_continuous(name ="Percentage of trees with Ctenophora-sister", breaks = seq(0,100,10), limits = c(0,110)) +
  scale_x_continuous(name = "Number of sites", breaks = seq(0,90000,15000), minor_breaks = seq(0,90000,5000)) +
  theme_bw()
# Number of constant sites
ggplot(data = plot_df, aes(x = proportion_constant_sites, y = value)) +
  geom_point() +
  scale_y_continuous(name ="Percentage of trees with Ctenophora-sister", breaks = seq(0,100,10), limits = c(0,110)) +
  scale_x_continuous(name = "Proportion of constant sites", breaks = seq(0.1,0.35,0.05), limits = c(0.1,0.35)) +
  theme_bw()
# Number of informative sites
ggplot(data = plot_df, aes(x = proportion_informative_sites, y = value)) +
  geom_point() +
  scale_y_continuous(name ="Percentage of trees with Ctenopora-sister", breaks = seq(0,100,10), limits = c(0,110)) +
  scale_x_continuous(name = "Proportion of informative sites", breaks = seq(0.5,0.8,0.05)) +
  theme_bw()

### Plot 4: Percent of trees with Ctenophora sister against length of Ctenophora branch lengths ###
plot_df <- melt(al_df, 
                id.vars = c("ID", "best_model", "num_taxa", "num_sites", "ctenophora_branch_length"),
                measure.vars = c("percent_CTEN_sister") )
ggplot(data = plot_df, aes(x = ctenophora_branch_length, y = value)) +
  geom_point() +
  scale_y_continuous(name ="Percentage of trees with Ctenophora-sister", breaks = seq(0,100,10), limits = c(0,100)) +
  scale_x_continuous(name = "Length of branch to Ctenophora clade", breaks = seq(0,0.80,0.1), minor_breaks = seq(0,0.80,0.05)) +
  theme_bw()

### Plot 5: Percent of trees with Ctenophora sister against length of Porifera branch lengths ###
plot_df <- melt(al_df, 
                id.vars = c("ID", "best_model", "num_taxa", "num_sites", "porifera_branch_length"),
                measure.vars = c("percent_CTEN_sister") )
# Number of sites
ggplot(data = plot_df, aes(x = porifera_branch_length, y = value)) +
  geom_point() +
  scale_y_continuous(name ="Percentage of trees with Ctenophora-sister", breaks = seq(0,100,10), limits = c(0,100)) +
  scale_x_continuous(name = "Length of branch to Porifera clade", breaks = seq(0,0.05,0.01), minor_breaks = seq(0,0.05,0.005)) +
  theme_bw()

### Plot 5: Percent of trees with Ctenophora sister length of Porifera branch lengths ###
plot_df <- melt(al_df, 
                id.vars = c("ID", "best_model", "num_taxa", "num_sites", "porifera_branch_length"),
                measure.vars = c("percent_CTEN_sister") )
# Number of sites
ggplot(data = plot_df, aes(x = porifera_branch_length, y = value)) +
  geom_point() +
  scale_y_continuous(name ="Percentage of trees with Ctenophora-sister", breaks = seq(0,100,10), limits = c(0,100)) +
  scale_x_continuous(name = "Length of branch to Porifera clade", breaks = seq(0,0.05,0.01), minor_breaks = seq(0,0.05,0.005)) +
  theme_bw()



#### 4. Plot tree weights from MAST model ####
# List all output files
all_files <- list.files(results_dir, recursive = TRUE)
mast_df_file <- paste0(results_dir, grep("summary_MAST_treeWeight_results", all_files, value = TRUE))
mast_df <- read.csv(mast_df_file, header = TRUE)
# Convert MAST output to long format
mast_long <- melt(mast_df,
                  id.vars = c("dataset", "matrix_name", "model_class", "model_code", "mast_branch_type", "minimum_branch_length", "number_hypothesis_trees", "year"),
                  measure.vars = c("tree_1_tree_weight", "tree_2_tree_weight"))
mast_long$var_label <- factor(mast_long$variable,
                              levels = c("tree_1_tree_weight", "tree_2_tree_weight"),
                              labels = c("Ctenophora-sister", "Porifera-sister"),
                              ordered = TRUE)
mast_long$dataset_label <- factor(mast_long$matrix_name,
                                  levels = c( "Dunn2008_FixedNames", "Philippe_etal_superalignment_FixedNames", "Pick2010",
                                              "UPDUNN_MB_FixedNames", "nonribosomal_9187_smatrix", "ribosomal_14615_smatrix",
                                              "REA_EST_includingXenoturbella", "ED3d", "Best108", "Chang_AA", 
                                              "Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned", "Metazoa_Choano_RCFV_strict",                             
                                              "Tplx_phylo_d1", "nonbilateria_MARE_BMGE"),
                                  labels = c("Dunn 2008", "Philippe 2009", "Pick 2010", "Philippe 2011", "Nosenko 2013\nnon-ribosomal", 
                                             "Nosenko 2013\nribosomal", "Ryan 2013", "Moroz 2014", "Borowiec 2015", "Chang 2015", 
                                             "Whelan 2015", "Whelan 2017", "Laumer 2018", "Laumer 2019" ),
                                  ordered = TRUE)
# Plot with boxplot for each dataset
bp <- ggplot(mast_long, aes(x = var_label, y = value, fill = var_label)) +
  geom_boxplot() +
  facet_wrap(~dataset_label) +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = "Tree weight", limits = c(0,1), breaks = seq(0,1,0.2), labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1)) +
  scale_fill_manual(name = "Hypothesis tree", labels = c("Ctenophora-sister", "Porifera-sister"), values = tree2_tonal, guide = "none") +
  labs(title = "MAST tree weights") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 25, margin = margin(t = 0, r = 15, b = 0, l = 10)),
        axis.text.x = element_text(size = 20, vjust = 0.5, hjust = 1, angle = 90, margin = margin(t = 10, r = 0, b = 10, l = 0)),  
        axis.text.y = element_text(size = 15),
        strip.text = element_text(size = 20),
        plot.title = element_text(size = 40, hjust = 0.5, margin = margin(t = 10, r = 0, b = 15, l = 0)) )
bp_file <- paste0(plot_dir, "MAST_tree_weights_2tree.")
ggsave(filename = paste0(bp_file, "png"), plot = bp, device = "png", width = 10, height = 14, units = "in")
ggsave(filename = paste0(bp_file, "pdf"), plot = bp, device = "pdf", width = 10, height = 14, units = "in")



#### 5. Plot expected likelihood weights from tree topology tests ####
# List all output files
all_files <- list.files(results_dir, recursive = TRUE)
au_df_file <- paste0(results_dir, grep("summary_au_test_results", all_files, value = TRUE))
au_df <- read.csv(au_df_file, header = TRUE)
# Convert MAST output to long format
au_long <- melt(au_df,
                id.vars = c("dataset", "matrix", "model_class", "best_model_code", "topology_test", "year"),
                measure.vars = c("tree_1", "tree_2"))
au_long$var_label <- factor(au_long$variable,
                            levels = c("tree_1", "tree_2"),
                            labels = c("Ctenophora-sister", "Porifera-sister"),
                            ordered = TRUE)
au_long$dataset_label <- factor(au_long$matrix,
                                levels = c( "Dunn2008_FixedNames", "Philippe_etal_superalignment_FixedNames", "Pick2010",
                                            "UPDUNN_MB_FixedNames", "nonribosomal_9187_smatrix", "ribosomal_14615_smatrix",
                                            "REA_EST_includingXenoturbella", "ED3d", "Best108", "Chang_AA", 
                                            "Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned", "Metazoa_Choano_RCFV_strict",                             
                                            "Tplx_phylo_d1", "nonbilateria_MARE_BMGE"),
                                labels = c("Dunn 2008", "Philippe 2009", "Pick 2010", "Philippe 2011", "Nosenko 2013\nnon-ribosomal", 
                                           "Nosenko 2013\nribosomal", "Ryan 2013", "Moroz 2014", "Borowiec 2015", "Chang 2015", 
                                           "Whelan 2015", "Whelan 2017", "Laumer 2018", "Laumer 2019" ),
                                ordered = TRUE)
# Plot with boxplot for each dataset
bp <- ggplot(au_long, aes(x = var_label, y = value, fill = var_label)) +
  geom_hline(yintercept = 0.05, color = "darkgrey", linetype = "dashed") +
  geom_boxplot() +
  facet_wrap(~dataset_label) +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = "p-value", limits = c(0,1), breaks = seq(0,1,0.2), labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1)) +
  scale_fill_manual(name = "Hypothesis tree", labels = c("Ctenophora-sister", "Porifera-sister"), values = tree2_tonal, guide = "none") +
  labs(title = "AU Test") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 25, margin = margin(t = 0, r = 15, b = 0, l = 10)),
        axis.text.x = element_text(size = 20, vjust = 0.5, hjust = 1, angle = 90, margin = margin(t = 10, r = 0, b = 10, l = 0)),  
        axis.text.y = element_text(size = 15),
        strip.text = element_text(size = 20),
        plot.title = element_text(size = 40, hjust = 0.5, margin = margin(t = 10, r = 0, b = 15, l = 0)) )
bp_file <- paste0(plot_dir, "au_test_2tree.")
ggsave(filename = paste0(bp_file, "png"), plot = bp, device = "png", width = 10, height = 14, units = "in")
ggsave(filename = paste0(bp_file, "pdf"), plot = bp, device = "pdf", width = 10, height = 14, units = "in")



#### 6. Plot expected likelihood weights from tree topology tests ####
# List all output files
all_files <- list.files(results_dir, recursive = TRUE)
elw_df_file <- paste0(results_dir, grep("summary_elw_results", all_files, value = TRUE))
elw_df <- read.csv(elw_df_file, header = TRUE)
# Convert MAST output to long format
elw_long <- melt(elw_df,
                 id.vars = c("dataset", "matrix", "model_class", "best_model_code", "topology_test", "year"),
                 measure.vars = c("tree_1", "tree_2"))
elw_long$var_label <- factor(elw_long$variable,
                             levels = c("tree_1", "tree_2"),
                             labels = c("Ctenophora-sister", "Porifera-sister"),
                             ordered = TRUE)
elw_long$dataset_label <- factor(elw_long$matrix,
                                 levels = c( "Dunn2008_FixedNames", "Philippe_etal_superalignment_FixedNames", "Pick2010",
                                             "UPDUNN_MB_FixedNames", "nonribosomal_9187_smatrix", "ribosomal_14615_smatrix",
                                             "REA_EST_includingXenoturbella", "ED3d", "Best108", "Chang_AA", 
                                             "Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned", "Metazoa_Choano_RCFV_strict",                             
                                             "Tplx_phylo_d1", "nonbilateria_MARE_BMGE"),
                                 labels = c("Dunn 2008", "Philippe 2009", "Pick 2010", "Philippe 2011", "Nosenko 2013\nnon-ribosomal", 
                                            "Nosenko 2013\nribosomal", "Ryan 2013", "Moroz 2014", "Borowiec 2015", "Chang 2015", 
                                            "Whelan 2015", "Whelan 2017", "Laumer 2018", "Laumer 2019" ),
                                 ordered = TRUE)
# Plot with boxplot for each dataset
bp <- ggplot(elw_long, aes(x = var_label, y = value, fill = var_label)) +
  geom_boxplot() +
  facet_wrap(~dataset_label) +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = "Weight", limits = c(0,1), breaks = seq(0,1,0.2), labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1)) +
  scale_fill_manual(name = "Hypothesis tree", labels = c("Ctenophora-sister", "Porifera-sister"), values = tree2_tonal, guide = "none") +
  labs(title = "Expected likelihood weight") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 25, margin = margin(t = 0, r = 15, b = 0, l = 10)),
        axis.text.x = element_text(size = 20, vjust = 0.5, hjust = 1, angle = 90, margin = margin(t = 10, r = 0, b = 10, l = 0)),  
        axis.text.y = element_text(size = 15),
        strip.text = element_text(size = 20),
        plot.title = element_text(size = 40, hjust = 0.5, margin = margin(t = 10, r = 0, b = 15, l = 0)) )
bp_file <- paste0(plot_dir, "expected_likelihood_weights_2tree.")
ggsave(filename = paste0(bp_file, "png"), plot = bp, device = "png", width = 10, height = 14, units = "in")
ggsave(filename = paste0(bp_file, "pdf"), plot = bp, device = "pdf", width = 10, height = 14, units = "in")


