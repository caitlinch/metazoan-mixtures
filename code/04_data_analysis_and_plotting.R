## caitlinch/metazoan-mixtures/code/03_data_analysis_and_plotting.R
# This script performs data analysis and creates plots of methods and results
# Caitlin Cherryh 2023

## This script:
# 1. Plots figures for the introduction and methods section of the manuscript


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
  scale_y_continuous(name ="Percentage of trees with Ctenopora-sister", breaks = seq(0,100,10), limits = c(0,110)) +
  scale_x_continuous(name = "Number of sites", breaks = seq(0,90000,15000), minor_breaks = seq(0,90000,5000)) +
  theme_bw()
# Number of constant sites
ggplot(data = plot_df, aes(x = proportion_constant_sites, y = value)) +
  geom_point() +
  scale_y_continuous(name ="Percentage of trees with Ctenopora-sister", breaks = seq(0,100,10), limits = c(0,110)) +
  scale_x_continuous(name = "Proportion of constant sites", breaks = seq(0.1,0.35,0.05), limits = c(0.1,0.35)) +
  theme_bw()
# Number of informative sites
ggplot(data = plot_df, aes(x = proportion_informative_sites, y = value)) +
  geom_point() +
  scale_y_continuous(name ="Percentage of trees with Ctenopora-sister", breaks = seq(0,100,10), limits = c(0,110)) +
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
  scale_y_continuous(name ="Percentage of trees with Ctenopora-sister", breaks = seq(0,100,10), limits = c(0,110)) +
  scale_x_continuous(name = "Number of sites", breaks = seq(0,90000,15000), minor_breaks = seq(0,90000,5000)) +
  theme_bw()
# Number of constant sites
ggplot(data = plot_df, aes(x = proportion_constant_sites, y = value)) +
  geom_point() +
  scale_y_continuous(name ="Percentage of trees with Ctenopora-sister", breaks = seq(0,100,10), limits = c(0,110)) +
  scale_x_continuous(name = "Proportion of constant sites", breaks = seq(0.1,0.35,0.05), limits = c(0.1,0.35)) +
  theme_bw()
# Number of informative sites
ggplot(data = plot_df, aes(x = proportion_informative_sites, y = value)) +
  geom_point() +
  scale_y_continuous(name ="Percentage of trees with Ctenopora-sister", breaks = seq(0,100,10), limits = c(0,110)) +
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
  scale_y_continuous(name ="Percentage of trees with Ctenopora-sister", breaks = seq(0,100,10), limits = c(0,110)) +
  scale_x_continuous(name = "Number of sites", breaks = seq(0,90000,15000), minor_breaks = seq(0,90000,5000)) +
  theme_bw()
# Number of constant sites
ggplot(data = plot_df, aes(x = proportion_constant_sites, y = value)) +
  geom_point() +
  scale_y_continuous(name ="Percentage of trees with Ctenopora-sister", breaks = seq(0,100,10), limits = c(0,110)) +
  scale_x_continuous(name = "Proportion of constant sites", breaks = seq(0.1,0.35,0.05), limits = c(0.1,0.35)) +
  theme_bw()
# Number of informative sites
ggplot(data = plot_df, aes(x = proportion_informative_sites, y = value)) +
  geom_point() +
  scale_y_continuous(name ="Percentage of trees with Ctenopora-sister", breaks = seq(0,100,10), limits = c(0,110)) +
  scale_x_continuous(name = "Proportion of informative sites", breaks = seq(0.5,0.8,0.05)) +
  theme_bw()










