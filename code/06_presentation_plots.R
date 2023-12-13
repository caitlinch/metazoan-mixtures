## caitlinch/metazoan-mixtures/code/06_presentation_plots.R
# This script plots data nicely for presentation slides
# Caitlin Cherryh 2023

## This script:
# 1. Plots figures for presentations and conference talks


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

# List all the output files from the analysis
all_output_files <- list.files(results_dir)

# Colour palette
metazoan_palette <- c(Bilateria = "#CC79A7", Cnidaria = "#009E73", Ctenophora = "#56B4E9", Porifera = "#E69F00", Outgroup = "#999999")
plasma_palette <- c("#0d0887", "#ed7953", "#9c179e", "#a0da39", "#440154") # Order: CTEN, PORI, CTEN+PORI, CTEN&PARA PORI, PORI & PARA PORI


#### 3. Prepare data for plotting ####
# Open the MAST output
mast_output_path <- paste0(results_dir, grep("MAST_model_output", all_output_files, value = T))
mast_op_df <- read.csv(mast_output_path, stringsAsFactors = FALSE)
mast_2 <- mast_op_df[mast_op_df$hypothesis_tree_analysis == "2_trees", ]
mast_5 <- mast_op_df[mast_op_df$hypothesis_tree_analysis == "5_trees", ]
# Melt the MAST output into a long-format dataframe consisting of only tree weights
tw_df <- melt(mast_op_df, 
              id.vars = c("dataset", "matrix_name", "model_code", "model_type", "mast_branch_type", "minimum_branch_length", "number_hypothesis_trees"),
              measure.vars = c("tree_1_tree_weight", "tree_2_tree_weight", "tree_3_tree_weight", "tree_4_tree_weight", "tree_5_tree_weight"))
# Plot as a violin plot
p <- ggplot(tw_df, aes(x = variable, y = value, fill = variable)) + 
  geom_violin() +
  scale_y_continuous(name = "Tree weight", breaks = seq(0,1,0.2), labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(0,1)) + 
  scale_x_discrete(name = "Evolutionary hypothesis", labels = c("CTEN", "PORI", "CTEN+PORI", "CTEN\n Para. PORI", "PORI\n Para. PORI")) +
  scale_fill_manual(guide = "none") +
  labs(title = "Applying the MAST model", subtitle = "MAST applied with best AA model and best PMSF model for 14 empirical datasets") + 
  theme_light() +
  theme(axis.title.x = element_text(size = 40, margin = margin(t = 30, r = 0, b = 10, l = 0), color = "grey30"), 
        axis.title.y = element_text(size = 40, margin = margin(t = 0, r = 30, b = 0, l = 10), color = "grey30"), 
        axis.text.x = element_text(size = 25, vjust = 1, hjust = 1, angle = 45, color = "grey30"),
        axis.text.y = element_text(size = 30, color = "grey30"),
        plot.title = element_text(size = 50, margin = margin(t = 10, r = 0, b = 30, l = 0), hjust = 0.5, color = "grey30"),
        plot.subtitle = element_text(size = 20, margin = margin(t = 0, r = 0, b = 10, l = 0), color = "grey30"))
plot_path <- paste0(plot_dir, "MAST_5_tree_weights_violin")
ggsave(filename = paste0(plot_path, ".png"), plot = p, device = "png")
ggsave(filename = paste0(plot_path, ".pdf"), plot = p, device = "pdf")


#### 4. Plotting 2-tree mixture ####
# Melt the MAST output into a long-format dataframe consisting of only tree weights
tw2_df <- melt(mast_2, 
               id.vars = c("dataset", "matrix_name", "model_code", "model_class", "mast_branch_type", "minimum_branch_length", "number_hypothesis_trees"),
               measure.vars = c("tree_1_tree_weight", "tree_2_tree_weight"))
# Plot as a violin plot
p_2t <- ggplot(tw2_df, aes(x = variable, y = value, fill = variable)) + 
  geom_violin() +
  scale_y_continuous(name = "Tree weight", breaks = seq(0,1,0.2), labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(0,1)) + 
  scale_x_discrete(name = "Evolutionary hypothesis", labels = c("CTEN", "PORI")) +
  scale_fill_manual(name = "Evolutionary hypothesis", values = plasma_palette[1:2]) +
  labs(title = "2-tree MAST Model") + 
  guides(fill="none") +
  theme_light() +
  theme(axis.title.x = element_text(size = 40, margin = margin(t = 30, r = 0, b = 10, l = 0), color = "grey30"), 
        axis.title.y = element_text(size = 40, margin = margin(t = 0, r = 30, b = 0, l = 10), color = "grey30"), 
        axis.text.x = element_text(size = 25, vjust = 1, hjust = 1, angle = 45, color = "grey30"),
        axis.text.y = element_text(size = 30, color = "grey30"),
        plot.title = element_text(size = 50, margin = margin(t = 10, r = 0, b = 30, l = 0), hjust = 0.5, color = "grey30"),
        plot.subtitle = element_text(size = 20, margin = margin(t = 0, r = 0, b = 10, l = 0), color = "grey30"),
        axis.line = element_line(linewidth = 1, color = "grey70"), 
        axis.ticks = element_line(linewidth = 1, color = "grey70"),
        panel.grid.major = element_line(linewidth = 1),
        panel.grid.minor = element_line(linewidth = 0.8),
        panel.border = element_rect(linewidth = 2, fill = NA, color = "grey70"))
# Save plot
p_2t_path <- paste0(plot_dir, "MAST_2_tree_weights_violin")
ggsave(filename = paste0(p_2t_path, ".png"), plot = p_2t, device = "png")
ggsave(filename = paste0(p_2t_path, ".pdf"), plot = p_2t, device = "pdf")



#### 5. Plotting 5-tree mixture ####
# Melt the MAST output into a long-format dataframe consisting of only tree weights
tw5_df <- melt(mast_5, 
               id.vars = c("dataset", "matrix_name", "model_code", "model_class", "mast_branch_type", "minimum_branch_length", "number_hypothesis_trees"),
               measure.vars = c("tree_1_tree_weight", "tree_2_tree_weight", "tree_3_tree_weight", "tree_4_tree_weight", "tree_5_tree_weight"))
# Plot as a violin plot
p_5t <- ggplot(tw5_df, aes(x = variable, y = value, fill = variable)) + 
  geom_violin() +
  scale_y_continuous(name = "Tree weight", breaks = seq(0,1,0.2), labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(0,1)) + 
  scale_x_discrete(name = "Evolutionary hypothesis", labels = c("CTEN", "PORI", "CTEN+PORI", "CTEN\n Para. PORI", "PORI\n Para. PORI")) +
  scale_fill_manual(name = "Evolutionary hypothesis", values = plasma_palette) +
  labs(title = "5-tree MAST Model") + 
  guides(fill="none") +
  theme_light() +
  theme(axis.title.x = element_text(size = 40, margin = margin(t = 30, r = 0, b = 10, l = 0), color = "grey30"), 
        axis.title.y = element_text(size = 40, margin = margin(t = 0, r = 30, b = 0, l = 10), color = "grey30"), 
        axis.text.x = element_text(size = 25, vjust = 1, hjust = 1, angle = 45, color = "grey30"),
        axis.text.y = element_text(size = 30, color = "grey30"),
        plot.title = element_text(size = 50, margin = margin(t = 10, r = 0, b = 30, l = 0), hjust = 0.5, color = "grey30"),
        plot.subtitle = element_text(size = 20, margin = margin(t = 0, r = 0, b = 10, l = 0), color = "grey30"),
        axis.line = element_line(linewidth = 1, color = "grey70"), 
        axis.ticks = element_line(linewidth = 1, color = "grey70"),
        panel.grid.major = element_line(linewidth = 1),
        panel.grid.minor = element_line(linewidth = 0.8),
        panel.border = element_rect(linewidth = 2, fill = NA, color = "grey70"))
# Save plot
p_5t_path <- paste0(plot_dir, "MAST_5_tree_weights_violin")
ggsave(filename = paste0(p_5t_path, ".png"), plot = p_5t, device = "png")
ggsave(filename = paste0(p_5t_path, ".pdf"), plot = p_5t, device = "pdf")



