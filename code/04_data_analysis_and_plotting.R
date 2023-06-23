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
repo_dir <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"



#### 2. Prepare functions, variables and packages ####
# Open packages
library(ggplot2)
library(ggtree)
library(patchwork)

# Open function files
source(paste0(repo_dir,"code/func_plotting.R"))



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

### Plot number of sites/number of informative sites against proportion of trees with each topology ###








