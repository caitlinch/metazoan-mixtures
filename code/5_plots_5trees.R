## caitlinch/metazoan-mixtures/code/05_plots_5trees.R
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
library(readxl)

# Open function files
source(paste0(repo_dir,"code/func_plotting.R"))
source(paste0(repo_dir,"code/func_data_processing.R"))

# Specify colour palettes (some of these are unused)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
metazoan_palette <- c(A = "#CC79A7", B = "#009E73", C = "#56B4E9", D = "#E69F00", E = "#999999")
tree2_palette <- c("#F0F921FF", "#0D0887FF")
tree5_palette <- c("#e66101", "#fdb863", "#f7f7f7", "#b2abd2", "#5e3c99")
tree5_cividis <- c("#FDE725FF", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF")
tree2_cividis <- c(tree5_cividis[1], tree5_cividis[5])
tree2_tonal <- c("#bdd7e7", "#2171b5")



#### 3. Plot tree weights from MAST model ####
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
bp_file <- paste0(plot_dir, "MAST_tree_weights_5tree.")
ggsave(filename = paste0(bp_file, "png"), plot = bp, device = "png", width = 10, height = 14, units = "in")
ggsave(filename = paste0(bp_file, "pdf"), plot = bp, device = "pdf", width = 10, height = 14, units = "in")



#### 4. Plot expected likelihood weights from tree topology tests ####
# List all output files
all_files <- list.files(results_dir, recursive = TRUE)
au_df_file <- paste0(results_dir, grep("summary_au_test_results", all_files, value = TRUE))
au_df <- read.csv(au_df_file, header = TRUE)
# Convert AU test output to long format
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
bp_file <- paste0(plot_dir, "au_test_5tree.")
ggsave(filename = paste0(bp_file, "png"), plot = bp, device = "png", width = 10, height = 14, units = "in")
ggsave(filename = paste0(bp_file, "pdf"), plot = bp, device = "pdf", width = 10, height = 14, units = "in")



#### 5. Plot expected likelihood weights from tree topology tests ####
# List all output files
all_files <- list.files(results_dir, recursive = TRUE)
elw_df_file <- paste0(results_dir, grep("summary_elw_results", all_files, value = TRUE))
elw_df <- read.csv(elw_df_file, header = TRUE)
# Convert ELW output to long format
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
bp_file <- paste0(plot_dir, "expected_likelihood_weights_5tree.")
ggsave(filename = paste0(bp_file, "png"), plot = bp, device = "png", width = 10, height = 14, units = "in")
ggsave(filename = paste0(bp_file, "pdf"), plot = bp, device = "pdf", width = 10, height = 14, units = "in")





