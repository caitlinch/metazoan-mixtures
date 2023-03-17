# metazoan-mixtures/code/03_data_analysis_and_plotting.R
## This script performs data analysis and creates plots of methods and results
# Caitlin Cherryh, 2022

## This script:
# 1. Plots figures for the introduction and methods section of the manuscript



#### 1. Input parameters ####
## Specify parameters:
# output_dir          <- Directory for IQ-Tree output (trees and tree mixtures)
# results_dir         <- Directory for results and plots
# repo_dir            <- Location of caitlinch/metazoan-mixtures github repository

output_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/05_plotting/"
repo_dir <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"



#### 2. Prepare functions, variables and packages ####
# Open packages
library(ggplot2)
library(ggtree)
library(patchwork)

# Open function files
source(paste0(repo_dir,"code/func_plotting.R"))



#### 3. Plot figures for introduction and methods sections of manuscript ####
## Plotting the alternative phylogenetic hypotheses ##
# Open the trees
trees_file <- paste0(repo_dir, "trees/alternative_phylogenetic_hypotheses.nex")
trees <- read.tree(trees_file)

# Plot one tree at a time
hyp1_plot <- ggtree(trees[[1]]) + geom_tree() + 
  geom_tiplab(geom = "text", size = 10) +
  theme_tree2(plot.margin=margin(6, 145, 6, 6)) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"))

hyp2_plot <- ggtree(trees[[2]]) + geom_tree() + 
  geom_tiplab(geom = "text", size = 10) +
  theme_tree2(plot.margin=margin(6, 145, 6, 6)) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"))

hyp3_plot <- ggtree(trees[[3]]) + geom_tree() + 
  geom_tiplab(geom = "text", size = 10) +
  theme_tree2(plot.margin=margin(6, 145, 6, 6)) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"))

hyp4_plot <- ggtree(trees[4]) + geom_tree() +
  geom_tiplab(geom = "text", size = 8) +
  theme_tree2(plot.margin=margin(6, 70, 6, 6)) +
  coord_cartesian(clip = 'off') +
  geom_cladelab(node=12, label="Porifera", textcolor="grey50", barcolor = "grey50", offset=7, geom = "text", align=TRUE, fontsize = 7, hjust = -0.1) +
  geom_cladelab(node=14, label="Porifera", textcolor="grey50", barcolor = "grey50", offset=7, geom = "text", align=TRUE, fontsize = 7, hjust = -0.1) +
  theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"))


hyp5_plot <- ggtree(trees[[5]]) + geom_tree() + 
  geom_tiplab(geom = "text", size = 8) +
  theme_tree2(plot.margin=margin(6, 70, 6, 6)) +
  coord_cartesian(clip = 'off') +
  geom_cladelab(node=11, label="Porifera", textcolor="grey50", barcolor = "grey50", offset=7, geom = "text", align=TRUE, fontsize = 7, hjust = -0.1) +
  geom_cladelab(node=13, label="Porifera", textcolor="grey50", barcolor = "grey50", offset=7, geom = "text", align=TRUE, fontsize = 7, hjust = -0.1) +
  theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"))

# Combine the four hypothesis trees into one 
patchwork_hyps <- (hyp1_plot | hyp2_plot | hyp3_plot)/(hyp4_plot | hyp5_plot) + 
  plot_annotation(tag_levels = 'a', tag_suffix = ".") & theme(plot.tag = element_text(size = 20))

# Export hypothesis plots 
hypothesis_plot_file <- paste0(output_dir, "hypothesis_tree_example_plot.png")
png(filename = hypothesis_plot_file, width = 1090, height = 723, units = "px", pointsize = 12, bg = "white")
patchwork_hyps
dev.off()



