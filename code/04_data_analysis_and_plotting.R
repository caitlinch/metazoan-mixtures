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
## Plotting the alternative phylogenetic hypotheses
# Open the trees
trees_file <- paste0(repo_dir, "trees/alternative_phylogenetic_hypotheses.nex")
trees <- read.tree(trees_file)

# Plot one tree at a time
hyp1_plot <- ggplot(trees[[1]]) + geom_tree() + 
  geom_tiplab(geom = "text", size = 10) +
  theme_tree2(plot.margin=margin(6, 145, 6, 6)) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"))

hyp2_plot <- ggplot(trees[[2]]) + geom_tree() + 
  geom_tiplab(geom = "text", size = 10) +
  theme_tree2(plot.margin=margin(6, 145, 6, 6)) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"))

hyp3_plot <- ggplot(trees[[3]]) + geom_tree() + 
  geom_tiplab(geom = "text", size = 10) +
  theme_tree2(plot.margin=margin(6, 145, 6, 6)) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"))

hyp4_plot <- ""

ggplot(trees[[4]]) + geom_tree() + 
  geom_tiplab(geom = "text", size = 9) +
  theme_tree2(plot.margin=margin(6, 80, 6, 6)) +
  coord_cartesian(clip = 'off') +
  geom_hilight(node = 11, fill = "gold") +
  geom_hilight(node = 13, fill = "gold") +
  geom_cladelab(node=11, label="Porifera", color="blue", offset=8, geom = "text", align=TRUE, fontsize = 7) +
  geom_cladelab(node=13, label="Porifera", color="blue", offset=8, geom = "text", align=TRUE, fontsize = 7) +
  theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"))


hyp5_plot <- ggplot(trees[[5]]) + geom_tree() + 
  geom_tiplab(geom = "text", size = 10) +
  theme_tree2(plot.margin=margin(6, 250, 6, 6)) +
  coord_cartesian(clip = 'off') +
  geom_hilight(node = 14, fill = "gold") +
  geom_hilight(node = 12, fill = "gold") +
  theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"))









