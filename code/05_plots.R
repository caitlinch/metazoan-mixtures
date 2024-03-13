## caitlinch/metazoan-mixtures/code/05_plots.R
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

## Control parameters
control_parameters <- list(add.extra.color.palettes = FALSE,
                           plot.hypothesis.trees = FALSE,
                           plot.MAST = TRUE,
                           plot.AU.tests = TRUE,
                           plot.ELW = FALSE,
                           plot.ML.topologies = TRUE,
                           plot.Porifera.topologies = FALSE)



#### 2. Prepare functions, variables and packages ####
# Open packages
library(ape)
library(ggplot2)
library(ggtree)
library(patchwork)
library(reshape2)
library(readxl)
library(dplyr)

# Open function files
source(paste0(repo_dir,"code/func_plotting.R"))
source(paste0(repo_dir,"code/func_data_processing.R"))

# Specify colour palettes used within these plots
metazoan_palette <- c(A = "#CC79A7", B = "#009E73", C = "#56B4E9", D = "#E69F00", E = "#999999")
model_class_qual <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")
names(model_class_qual) <- c("PM", "PMSF", "Mixture", "Q")
# Notes for Viridis color palette usage (function = scale_color_viridis_d):
#   For ML tree topology (i.e. Sister to all other Metazoans) use option = "C"/"plasma"
#   For Porifera clade topology (i.e. Monophyletic/Paraphyletic) use option = "D"/"viridis"

# Extra colour palettes (unused)
if (control_parameters$add.extra.color.palettes == TRUE){
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  metazoan_clade_palette <- c(Bilateria = "#CC79A7", Cnidaria = "#009E73", Ctenophora = "#56B4E9", Porifera = "#E69F00", Outgroup = "#999999")
  tree2_palette <- c("#F0F921FF", "#0D0887FF")
  tree5_palette <- c("#e66101", "#fdb863", "#f7f7f7", "#b2abd2", "#5e3c99")
  tree5_cividis <- c("#FDE725FF", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF")
  tree2_cividis <- c(tree5_cividis[1], tree5_cividis[5])
  tree2_tonal <- c("#bdd7e7", "#2171b5")
  model3_tonal <- c("#980043", "#df65b0", "#d4b9da")
}

# List all files
all_files <- list.files(results_dir, recursive = TRUE)
# Remove any for the 5 tree model
all_files <- grep("5trees|5_trees", all_files, value = TRUE, invert = TRUE) 
# Remove any excel files into separate object
excel_files <- grep("xls|xlsx", all_files, value = TRUE) 
all_files <- grep("xls|xlsx", all_files, value = TRUE, invert = TRUE) 



#### 3. Plot figures for introduction and methods sections of manuscript ####
if (control_parameters$plot.hypothesis.trees == TRUE){
  ### Plotting the alternative phylogenetic hypotheses ###
  ## Open the trees
  trees_file <- paste0(repo_dir, "trees/alternative_phylogenetic_hypotheses.nex")
  trees <- read.tree(trees_file)
  image_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/99_RSB_HDR_conference_2023/pictures/"
  
  ## Create a dataframe for labelling
  metazoan_clade_labs <- data.frame("clade" = c("Bilateria", "Cnidaria", "Ctenophora", "Porifera", "Calcarea", "Demospongiae", "Hexactinellida", "Homoscleromorpha", "Outgroup"),
                                    "color" = c("A", "B", "C", rep("D", 5), "E"))
  metazoan_clade_labs$lab <- metazoan_clade_labs$clade
  metazoan_clade_labs$html_lab <- paste0("<b style='color:", metazoan_clade_labs$color, "'>", metazoan_clade_labs$clade, "</b>")
  
  ## Plot the hypothesis trees
  ## Black and White plots: Combine the five hypothesis trees into one plot
  patchwork_hyps <- (bw.monophyletic.clades.plot(trees[[1]]) | bw.monophyletic.clades.plot(trees[[2]]) | bw.monophyletic.clades.plot(trees[[3]]))/
    (bw.paraphyletic.clades.plot(trees[[4]], label_nodes = c(12, 14)) | bw.paraphyletic.clades.plot(trees[[5]], label_nodes = c(11, 13))) + 
    plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 20))
  # Export hypothesis plot as png
  hypothesis_plot_file <- paste0(plot_dir, "hypothesis_tree_example_plot.")
  png(filename = paste0(hypothesis_plot_file, "png"), width = 1090, height = 723, units = "px", pointsize = 12, bg = "white")
  patchwork_hyps
  dev.off()
  # Export hypothesis plot as pdf
  pdf(file = paste0(hypothesis_plot_file, "pdf"), width = 15, height = 8)
  patchwork_hyps
  dev.off()
  
  ## Colour plots: Combine the five hypothesis trees into one plot 
  # Create individual plots
  p1 = color.clades.plot(trees[[1]], tip_labels = metazoan_clade_labs, color_palette = metazoan_palette, xlimits = c(0,5.5))
  p2 = color.clades.plot(trees[[2]], tip_labels = metazoan_clade_labs, color_palette = metazoan_palette, xlimits = c(0,5.5))
  p3 = color.clades.plot(trees[[3]], tip_labels = metazoan_clade_labs, color_palette = metazoan_palette, xlimits = c(0,4))
  p4 = color.clades.plot(trees[[4]], tip_labels = metazoan_clade_labs, color_palette = metazoan_palette, xlimits = c(0,9))
  p5 = color.clades.plot(trees[[5]], tip_labels = metazoan_clade_labs, color_palette = metazoan_palette, xlimits = c(0,9))
  # Collate plots using patchwork
  patchwork_hyps_color <- wrap_plots(p1, p2) / wrap_plots(p3, p4) / wrap_plots(p5, plot_spacer()) +
    plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 30)) 
  # Export hypothesis plot as png
  hypothesis_plot_file <- paste0(plot_dir, "hypothesis_tree_example_plot_color.")
  png(filename = paste0(hypothesis_plot_file, "png"), width = 1200, height = 900, units = "px", pointsize = 12, bg = "white")
  patchwork_hyps_color
  dev.off()
  # Export hypothesis plot as pdf
  pdf(file = paste0(hypothesis_plot_file, "pdf"), width = 17, height = 12)
  patchwork_hyps_color
  dev.off()
}



#### 3. Plot tree weights from MAST model ####
if (control_parameters$plot.MAST == TRUE){
  # Identify MAST csv file
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
  mast_long$model_class <- factor(mast_long$model_class,
                                  levels = c("CXX", "PMSF", "Other", "Single"),
                                  labels = c("PM", "PMSF", "Mixture", "Q"),
                                  ordered = T)
  # Plot with lines for each dataset/model class
  bp <- ggplot(mast_long, aes(x = var_label, y = value, color = model_class, group = model_class)) +
    geom_point(size = 3, alpha = 0.6) +
    geom_line(alpha = 0.6) +
    facet_wrap(~dataset_label) +
    scale_x_discrete(name = NULL) +
    scale_y_continuous(name = "Tree weight", limits = c(0,1), breaks = seq(0,1,0.2), labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1)) +
    labs(title = "MAST tree weights") +
    scale_color_manual(name = "Model class", values = model_class_qual) +
    theme_bw() +
    theme(axis.title.y = element_text(size = 25, margin = margin(t = 0, r = 15, b = 0, l = 10)),
          axis.text.x = element_text(size = 17, vjust = 0.5, hjust = 1, angle = 90, margin = margin(t = 10, r = 0, b = 10, l = 0)),  
          axis.text.y = element_text(size = 17),
          strip.text = element_text(size = 20),
          plot.title = element_text(size = 40, hjust = 0.5, margin = margin(t = 10, r = 0, b = 15, l = 0)),
          legend.title = element_text(size = 28),
          legend.text = element_text(size = 25) )
  bp_file <- paste0(plot_dir, "mainfigure_MAST_tree_weights_2tree.")
  ggsave(filename = paste0(bp_file, "png"), plot = bp, device = "png", width = 12, height = 14, units = "in")
  ggsave(filename = paste0(bp_file, "pdf"), plot = bp, device = "pdf", width = 12, height = 14, units = "in")
 }



#### 4. Plot AU test results from tree topology tests ####
if (control_parameters$plot.AU.test == TRUE){
  # Identify AU test results
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
  au_long$model_class <- factor(au_long$best_model_code,
                                levels = c("LG_C60", "LG_C20", "C60", "C20",
                                           "PMSF_C60", "PMSF_LG_C60",
                                           "UL3", "LG4M",
                                           "GTR20_no_I"),
                                labels = c("PM", "PM", "PM", "PM",
                                           "PMSF", "PMSF",
                                           "Mixture", "Mixture", 
                                           "Q"),
                                ordered = T)
  # Plot with lines for each dataset/model class
  bp <- ggplot(au_long, aes(x = var_label, y = value, color = model_class, group = model_class)) +
    geom_hline(yintercept = 0.05, color = "darkgrey", linetype = "dashed") +
    geom_point(size = 3, alpha = 0.6) +
    geom_line(alpha = 0.6) +
    facet_wrap(~dataset_label) +
    scale_x_discrete(name = NULL) +
    scale_y_continuous(name = "p-value", limits = c(0,1), breaks = seq(0,1,0.2), labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1)) +
    scale_color_manual(name = "Model class", values = model_class_qual) +
    labs(title = "AU Test") +
    theme_bw() +
    theme(axis.title.y = element_text(size = 25, margin = margin(t = 0, r = 15, b = 0, l = 10)),
          axis.text.x = element_text(size = 17, vjust = 0.5, hjust = 1, angle = 90, margin = margin(t = 10, r = 0, b = 10, l = 0)),  
          axis.text.y = element_text(size = 17),
          strip.text = element_text(size = 20),
          plot.title = element_text(size = 40, hjust = 0.5, margin = margin(t = 10, r = 0, b = 15, l = 0)),
          legend.title = element_text(size =  28),
          legend.text = element_text(size = 25) )
  bp_file <- paste0(plot_dir, "mainfigure_au_test_2tree.")
  ggsave(filename = paste0(bp_file, "png"), plot = bp, device = "png", width = 12, height = 14, units = "in")
  ggsave(filename = paste0(bp_file, "pdf"), plot = bp, device = "pdf", width = 12, height = 14, units = "in")
 }



#### 5. Plot expected likelihood weights from tree topology tests ####
if (control_parameters$plot.ELW == TRUE){
  # Identify ELW results
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
  elw_long$model_class <- factor(elw_long$model_class,
                                 levels = c("LG_C60", "LG_C20", "C60", "C20",
                                            "PMSF_C60", "PMSF_LG_C60",
                                            "UL3", "LG4M",
                                            "GTR20_no_I"),
                                 labels = c("PM", "PM", "PM", "PM",
                                            "PMSF", "PMSF",
                                            "Mixture", "Mixture", 
                                            "Q"),
                                 ordered = T)
  # Plot with lines for each dataset/model class
  bp <- ggplot(elw_long, aes(x = var_label, y = value, color = model_class, group = model_class)) +
    geom_point(size = 3, alpha = 0.6) + 
    geom_line(alpha = 0.6) +
    facet_wrap(~dataset_label) +
    scale_x_discrete(name = NULL) +
    scale_y_continuous(name = "Weight", limits = c(0,1), breaks = seq(0,1,0.2), labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1)) +
    scale_color_manual(name = "Model class", values = model_class_qual) +
    labs(title = "Expected likelihood weight") +
    theme_bw() +
    theme(axis.title.y = element_text(size = 25, margin = margin(t = 0, r = 15, b = 0, l = 10)),
          axis.text.x = element_text(size = 15, vjust = 0.5, hjust = 1, angle = 90, margin = margin(t = 10, r = 0, b = 10, l = 0)),  
          axis.text.y = element_text(size = 15),
          strip.text = element_text(size = 20),
          plot.title = element_text(size = 40, hjust = 0.5, margin = margin(t = 10, r = 0, b = 15, l = 0)),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15) )
  bp_file <- paste0(plot_dir, "mainfigure_expected_likelihood_weights_2tree.")
  ggsave(filename = paste0(bp_file, "png"), plot = bp, device = "png", width = 12, height = 14, units = "in")
  ggsave(filename = paste0(bp_file, "pdf"), plot = bp, device = "pdf", width = 12, height = 14, units = "in")
}



#### 6. Summarise maximum likelihood topology tree results ####
if (control_parameters$plot.ML.topologies == TRUE | control_parameters$plot.Porifera.topologies == TRUE){
  ## Plot tree topology
  # Open dataframe
  topo_df_file <- grep("xls", grep("ML_tree_topology_ManualCheck", excel_files, value = TRUE), value = TRUE)
  topo_df_file <- grep("Summary", topo_df_file, value = TRUE, invert = TRUE)
  topo_df <- as.data.frame(read_excel(path = paste0(results_dir, "/", topo_df_file), sheet = "Topology"))
  # Convert topology output to long format
  topo_long <- melt(topo_df,
                    id.vars = c("dataset", "matrix_name", "model_code", "PORI_topology"),
                    measure.vars = c("sister_group"))
  topo_long$dataset_label <- factor(topo_long$matrix_name,
                                    levels = c("Dunn2008_FixedNames", "Philippe_etal_superalignment_FixedNames", "Pick2010",
                                               "UPDUNN_MB_FixedNames", "nonribosomal_9187_smatrix", "ribosomal_14615_smatrix",
                                               "REA_EST_includingXenoturbella", "ED3d", "Best108", "Chang_AA",  "Dataset10",
                                               "Metazoa_Choano_RCFV_strict", "Tplx_phylo_d1", "nonbilateria_MARE_BMGE"),
                                    labels = c("Dunn 2008", "Philippe 2009", "Pick 2010", "Philippe 2011", "Nosenko 2013\nnon-ribosomal", 
                                               "Nosenko 2013\nribosomal", "Ryan 2013", "Moroz 2014", "Borowiec 2015", "Chang 2015", 
                                               "Whelan 2015", "Whelan 2017", "Laumer 2018", "Laumer 2019" ),
                                    ordered = TRUE)
  topo_long$dataset_label_singleLine <- factor(topo_long$matrix_name,
                                               levels = c("Dunn2008_FixedNames", "Philippe_etal_superalignment_FixedNames", "Pick2010",
                                                          "UPDUNN_MB_FixedNames", "nonribosomal_9187_smatrix", "ribosomal_14615_smatrix",
                                                          "REA_EST_includingXenoturbella", "ED3d", "Best108", "Chang_AA",  "Dataset10",
                                                          "Metazoa_Choano_RCFV_strict", "Tplx_phylo_d1", "nonbilateria_MARE_BMGE"),
                                               labels = c("Dunn 2008", "Philippe 2009", "Pick 2010", "Philippe 2011", "Nosenko 2013 non-ribosomal", 
                                                          "Nosenko 2013 ribosomal", "Ryan 2013", "Moroz 2014", "Borowiec 2015", "Chang 2015", 
                                                          "Whelan 2015", "Whelan 2017", "Laumer 2018", "Laumer 2019" ),
                                               ordered = TRUE)
  topo_long$PORI_topology <- factor(topo_long$PORI_topology,
                                    levels = c("One taxon", "Monophyletic", "Paraphyletic"),
                                    labels = c("One taxon", "Monophyletic", "Paraphyletic"),
                                    ordered = TRUE)
  
  ## Plot number of models with each topology
  # Plot with barchart for each dataset - one bar per dataset
  bc <- ggplot(topo_long, aes(x = dataset_label_singleLine, fill = value)) +
    geom_bar() +
    labs(title = "ML tree topology") +
    scale_x_discrete(name = NULL) +
    scale_y_continuous(name = "Number of models", limits = c(0,26), breaks = seq(0,30,4), labels = seq(0,30,4), minor_breaks = seq(0,30,2)) +
    scale_fill_viridis_d(name = "Sister to other\nMetazoan clades", option = "C") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 25, margin = margin(t = 0, r = 15, b = 0, l = 10)),
          axis.text.x = element_text(size = 15, hjust = 1, vjust = 0.5, angle = 90),
          axis.text.y = element_text(size = 15),
          strip.text = element_text(size = 20),
          plot.title = element_text(size = 30, hjust = 0.5, margin = margin(t = 10, r = 0, b = 15, l = 0)),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15))
  bc_file <- paste0(plot_dir, "ML_topology_results_singleBar.")
  ggsave(filename = paste0(bc_file, "png"), plot = bc, device = "png")
  ggsave(filename = paste0(bc_file, "pdf"), plot = bc, device = "pdf")
  
  ## Plot number of each models with each Porifera clade topology
  # Convert topology output to long format
  bc2 <- ggplot(topo_long, aes(x = dataset_label_singleLine, fill = PORI_topology)) +
    geom_bar() +
    labs(title = "Porifera topology") +
    scale_x_discrete(name = NULL) +
    scale_y_continuous(name = "Number of models", limits = c(0,26), breaks = seq(0,30,4), labels = seq(0,30,4), minor_breaks = seq(0,30,2)) +
    scale_fill_viridis_d(name = "Porifera clade\ntopology", option = "D") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 25, margin = margin(t = 0, r = 15, b = 0, l = 10)),
          axis.text.x = element_text(size = 15, hjust = 1, vjust = 0.5, angle = 90),
          axis.text.y = element_text(size = 15),
          strip.text = element_text(size = 20),
          plot.title = element_text(size = 30, hjust = 0.5, margin = margin(t = 10, r = 0, b = 15, l = 0)),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15))
  bc2_file <- paste0(plot_dir, "Porifera_topology_results_singleBar.")
  ggsave(filename = paste0(bc2_file, "png"), plot = bc2, device = "png")
  ggsave(filename = paste0(bc2_file, "pdf"), plot = bc2, device = "pdf")
  
  ## Combined plots of number of models and number of Porifera topologies
  # Plot ML tree topology
  bc <- ggplot(topo_long, aes(x = dataset_label_singleLine, fill = value)) +
    geom_bar() +
    scale_x_discrete(name = NULL) +
    scale_y_continuous(name = "Number of models", limits = c(0,26), breaks = seq(0,30,4), labels = seq(0,30,4), minor_breaks = seq(0,30,2)) +
    scale_fill_viridis_d(name = "Sister to other\nMetazoan clades", option = "C") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.text.x = element_text(size = 12, hjust = 1, vjust = 0.5, angle = 90),
          axis.text.y = element_text(size = 12),
          strip.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12))
  # Plot Porifera clade topology
  bc2 <- ggplot(topo_long, aes(x = dataset_label_singleLine, fill = PORI_topology)) +
    geom_bar() +
    scale_x_discrete(name = NULL) +
    scale_y_continuous(name = "Number of models", limits = c(0,26), breaks = seq(0,30,4), labels = seq(0,30,4), minor_breaks = seq(0,30,2)) +
    scale_fill_viridis_d(name = "Porifera clade\ntopology", option = "D") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.text.x = element_text(size = 12, hjust = 1, vjust = 0.5, angle = 90),
          axis.text.y = element_text(size = 12),
          strip.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12))
  # Combine into a quilt
  quilt <- (bc / bc2) + plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 20))
  # Save and output quilt
  quilt_file <- paste0(plot_dir, "mainfigure_combined_topology_results_singleBar.")
  ggsave(filename = paste0(quilt_file, "png"), plot = quilt, device = "png", units = "in", width = 8, height = 10)
  ggsave(filename = paste0(quilt_file, "pdf"), plot = quilt, device = "pdf", units = "in", width = 8, height = 10)
}



#### 7. Summarise maximum likelihood topology tree results by model class ####
if (control_parameters$plot.ML.topologies == TRUE | control_parameters$plot.Porifera.topologies == TRUE){
  ## Plot tree topology
  # Open dataframe
  topo_df_file <- grep("xls", grep("ML_tree_topology_ManualCheck", excel_files, value = TRUE), value = TRUE)
  topo_df_file <- grep("Summary", topo_df_file, value = TRUE, invert = TRUE)
  topo_df <- as.data.frame(read_excel(path = paste0(results_dir, "/", topo_df_file), sheet = "Topology"))
  # Remove ModelFinder row
  topo_df <- topo_df[topo_df$model_code != "ModelFinder",]
  # Convert topology output to long format
  topo_long <- melt(topo_df,
                    id.vars = c("dataset", "matrix_name", "model_code", "PORI_topology"),
                    measure.vars = c("sister_group"))
  topo_long$dataset_label <- factor(topo_long$matrix_name,
                                    levels = c("Dunn2008_FixedNames", "Philippe_etal_superalignment_FixedNames", "Pick2010",
                                               "UPDUNN_MB_FixedNames", "nonribosomal_9187_smatrix", "ribosomal_14615_smatrix",
                                               "REA_EST_includingXenoturbella", "ED3d", "Best108", "Chang_AA",  "Dataset10",
                                               "Metazoa_Choano_RCFV_strict", "Tplx_phylo_d1", "nonbilateria_MARE_BMGE"),
                                    labels = c("Dunn 2008", "Philippe 2009", "Pick 2010", "Philippe 2011", "Nosenko 2013\nnon-ribosomal", 
                                               "Nosenko 2013\nribosomal", "Ryan 2013", "Moroz 2014", "Borowiec 2015", "Chang 2015", 
                                               "Whelan 2015", "Whelan 2017", "Laumer 2018", "Laumer 2019" ),
                                    ordered = TRUE)
  topo_long$dataset_label_singleLine <- factor(topo_long$matrix_name,
                                               levels = c("Dunn2008_FixedNames", "Philippe_etal_superalignment_FixedNames", "Pick2010",
                                                          "UPDUNN_MB_FixedNames", "nonribosomal_9187_smatrix", "ribosomal_14615_smatrix",
                                                          "REA_EST_includingXenoturbella", "ED3d", "Best108", "Chang_AA",  "Dataset10",
                                                          "Metazoa_Choano_RCFV_strict", "Tplx_phylo_d1", "nonbilateria_MARE_BMGE"),
                                               labels = c("Dunn 2008", "Philippe 2009", "Pick 2010", "Philippe 2011", "Nosenko 2013 non-ribosomal", 
                                                          "Nosenko 2013 ribosomal", "Ryan 2013", "Moroz 2014", "Borowiec 2015", "Chang 2015", 
                                                          "Whelan 2015", "Whelan 2017", "Laumer 2018", "Laumer 2019" ),
                                               ordered = TRUE)
  topo_long$PORI_topology <- factor(topo_long$PORI_topology,
                                    levels = c("One taxon", "Monophyletic", "Paraphyletic"),
                                    labels = c("One taxon", "Monophyletic", "Paraphyletic"),
                                    ordered = TRUE)
  # Add model class to the topology_df
  #   Q is for single matrix amino-acid exchange rate matrices
  #   Mixture is is for protein mixture models from IQ-Tree
  #   ModelFinder is the results from ModelFinder in IQ-Tree - could be single matrix or protein mixture model
  #   PM is for "empirical profile mixture models" i.e. C20 and C60 models
  #   PMSF is for "posterior site mean frequency" models in IQ-Tree2 i.e. PMSF+C20 and PMSF+C60 models
  topo_long$model_class <- factor(topo_long$model_code,
                                  levels = c("GTR20", "JTT", "JTTDCMut", "LG", "mtZOA", "PMB", "Poisson", "rtREV", "WAG",
                                             "CF4", "EHO", "EX_EHO", "EX2", "EX3", "LG4M", "UL2", "UL3",
                                             "PMSF_C20", "PMSF_C60", "PMSF_LG_C20", "PMSF_LG_C60",
                                             "C20", "C60", "LG_C20", "LG_C60"),
                                  labels = c(rep("Q", 9), rep("Mixture", 8), rep("PMSF", 4), rep("PM", 4)),
                                  ordered = TRUE)
  
  ## Plot number of models with each topology - facet per dataset, one bar per model_class
  # Plot with barchart for each dataset - one bar per dataset
  bc <- ggplot(topo_long, aes(x = model_class, fill = value)) +
    geom_bar() +
    facet_wrap(~dataset_label, ncol = 3) +
    labs(title = "Tree topology") +
    scale_x_discrete(name = NULL) +
    scale_y_continuous(name = "Number of models", limits = c(0,9), breaks = seq(0,9,3), labels = seq(0,9,3), minor_breaks = seq(0,10,1)) +
    scale_fill_viridis_d(name = "Sister to other\nMetazoan clades", option = "C") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 15, b = 0, l = 10)),
          axis.text.x = element_text(size = 14, hjust = 1, vjust = 0.5, angle = 90),
          axis.text.y = element_text(size = 14),
          strip.text = element_text(size = 15),
          plot.title = element_text(size = 30, hjust = 0.5, margin = margin(t = 10, r = 0, b = 15, l = 0)),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 14))
  bc_file <- paste0(plot_dir, "mainfigure_ML_topology_results_singleBar_datasetFacet.")
  ggsave(filename = paste0(bc_file, "png"), plot = bc, device = "png", height = 10, width = 8, units = "in")
  ggsave(filename = paste0(bc_file, "pdf"), plot = bc, device = "pdf", height = 10, width = 8, units = "in")
  
  ## Plot number of each models with each Porifera clade topology - facet per dataset, one bar per model_class
  # Convert topology output to long format
  bc2 <- ggplot(topo_long, aes(x = model_class, fill = PORI_topology)) +
    geom_bar() +
    facet_wrap(~dataset_label, ncol = 3) +
    labs(title = "Porifera clade topology") +
    scale_x_discrete(name = NULL) +
    scale_y_continuous(name = "Number of models", limits = c(0,9), breaks = seq(0,9,3), labels = seq(0,9,3), minor_breaks = seq(0,10,1)) +
    scale_fill_viridis_d(name = "Porifera clade\ntopology", option = "D") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 15, b = 0, l = 10)),
          axis.text.x = element_text(size = 14, hjust = 1, vjust = 0.5, angle = 90),
          axis.text.y = element_text(size = 14),
          strip.text = element_text(size = 15),
          plot.title = element_text(size = 30, hjust = 0.5, margin = margin(t = 10, r = 0, b = 15, l = 0)),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 14))
  bc2_file <- paste0(plot_dir, "mainfigure_Porifera_topology_results_singleBar_datasetFacet.")
  ggsave(filename = paste0(bc2_file, "png"), plot = bc2, device = "png", height = 10, width = 8, units = "in")
  ggsave(filename = paste0(bc2_file, "pdf"), plot = bc2, device = "pdf", height = 10, width = 8, units = "in")
}




