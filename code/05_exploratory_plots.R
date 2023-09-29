## caitlinch/metazoan-mixtures/code/05_exploratory_plots.R
# This script performs data analysis and creates plots of methods and results
# Caitlin Cherryh 2023

## This script:
# 1. Plots exploratory figures with the results from the maximum likelihood tree analysis


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

# Specify colour palettes (some of these are unused)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
metazoan_palette <- c(A = "#CC79A7", B = "#009E73", C = "#56B4E9", D = "#E69F00", E = "#999999")
metazoan_clade_palette <- c(Bilateria = "#CC79A7", Cnidaria = "#009E73", Ctenophora = "#56B4E9", Porifera = "#E69F00", Outgroup = "#999999")
tree2_palette <- c("#F0F921FF", "#0D0887FF")
tree5_palette <- c("#e66101", "#fdb863", "#f7f7f7", "#b2abd2", "#5e3c99")
tree5_cividis <- c("#FDE725FF", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF")
tree2_cividis <- c(tree5_cividis[1], tree5_cividis[5])
tree2_tonal <- c("#bdd7e7", "#2171b5")

# List all files
all_files <- list.files(results_dir, recursive = TRUE)
# Remove any for the 5 tree model
all_files <- grep("5trees", all_files, value = TRUE, invert = TRUE) 



#### 3. Exploratory plots ####
# Open the summary alignment details and the summary topology results
al_df_file_path <- paste0(results_dir, "04_02_exploratory_data_analysis_branch_lengths.csv")
if (file.exists(al_df_file_path) == TRUE){
  al_df <- read.csv(al_df_file_path, stringsAsFactors = FALSE)
} else {
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
  al_df$dataset_label <- c("Dunn2008", "Philippe2009", "Pick2010", "Philippe2011", "Nosenko2013 non-ribo",
                           "Nosenko2013 ribo", "Ryan2013", "Moroz2014", "Borowiec2015", "Chang2015",
                           "Whelan2015", "Whelan2017", "Laumer2018", "Laumer2019")
  # Replicate each row of the al_df three times
  al_df <- al_df %>% slice(rep(1:n(), each = 3))
  # Open the best_model df to copy the combinations of alignment and model
  best_model_df <- read.csv(paste0(results_dir, grep("MAST_parameters.csv", all_files, value = T)), stringsAsFactors = FALSE)
  best_model_df$matrix_name[which(best_model_df$dataset == "Whelan2015")] <- "Dataset10"
  # Order al_df to match the best_model df
  best_model_df <- best_model_df[order(best_model_df$dataset, best_model_df$matrix_name, best_model_df$model_class), ]
  al_df <- al_df[order(al_df$dataset, al_df$matrix_name), ]
  # Add a column for best model, model code and model class
  al_df$model_class <- best_model_df$model_class
  al_df$model_code <- best_model_df$model_code
  al_df$best_model <- gsub("'", "", best_model_df$best_model)
  al_df$ID <- paste0(al_df$dataset, ".", al_df$matrix_name, ".", al_df$model_code)
  # Extract branch a and branch b lengths
  bl_df <- as.data.frame(do.call(rbind, (lapply(1:nrow(al_df), extract.branch.length.wrapper, alignment_df = al_df, tree_directory = tree_dir))))
  al_df <- cbind(al_df, bl_df)
  # Reorder al_df 
  al_df <- al_df[, c("dataset", "dataset_label", "matrix_name", "model_class", "model_code", "best_model", "ID", "sequence_format", 
                     "num_taxa", "num_sites", "number_constant_sites", "proportion_constant_sites", "number_invariant_sites",
                     "proportion_invariant_sites", "number_informative_sites", "proportion_informative_sites", "percent_CTEN_sister",
                     "percent_PORI_sister", "percent_CTEN_PORI_sister", "percent_Radiata_sister", "percent_PORI_monophyletic",
                     "percent_PORI_paraphyletic", "percent_PORI_one_taxon", "percent_CTEN_CNID_monophyletic", "percent_CTEN_CNID_not_monophyletic",
                     "ctenophora_clade_branch_length", "ctenophora_clade_depth", "porifera_clade_branch_length", "porifera_clade_depth")]
  # Save this exploratory analysis dataframe
  write.csv(al_df, file = al_df_file_path, row.names = FALSE)
}

### Plot number of sites/number of informative sites against proportion of trees with each topology ###
### Plot 1: Percent of trees with Ctenophora sister against number of sites ###
# Create labeller function
var_labs <- c("Number of sites", "Proportion of\nconstant sites", "Proportion of\ninvariant sites", "Proportion of\ninformative sites")
names(var_labs) = c("num_sites", "proportion_constant_sites", "proportion_invariant_sites", "proportion_informative_sites")
# Create plot
plot_df <- melt(al_df, 
                id.vars = c("ID", "best_model", "num_taxa", "num_sites", "percent_CTEN_sister"),
                measure.vars = c("num_sites", "proportion_invariant_sites", "proportion_informative_sites") )
p <- ggplot(data = plot_df, aes(x = value, y = percent_CTEN_sister)) +
  facet_grid(.~variable, scales = "free_x", labeller = labeller(variable = var_labs)) +
  geom_point(size = 3, color = "darkblue", alpha = 0.4) +
  scale_y_continuous(name ="Percentage of trees with Ctenophora-sister", breaks = seq(0,100,10), limits = c(0,110)) +
  scale_x_continuous(name = "Number of sites / Proportion of sites") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 20, margin = margin(t = 15, r = 0, b = 10, l = 0)),
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 15, b = 0, l = 10)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        strip.text = element_text(size = 20))
p_path <- paste0(plot_dir, "exploratory_NumSites_Ctenophora-sister.")
ggsave(filename = paste0(p_path, "png"), plot = p, device = "png")
ggsave(filename = paste0(p_path, "pdf"), plot = p, device = "pdf")

### Plot 2: Percent of trees with monophyletic Porifer against number of sites ###
plot_df <- melt(al_df, 
                id.vars = c("ID", "best_model", "num_taxa", "num_sites", "percent_PORI_sister"),
                measure.vars = c("num_sites", "proportion_invariant_sites", "proportion_informative_sites") )
p <- ggplot(data = plot_df, aes(x = value, y = percent_PORI_sister)) +
  facet_grid(.~variable, scales = "free_x", labeller = labeller(variable = var_labs)) +
  geom_point(size = 3, color = "darkblue", alpha = 0.4) +
  scale_y_continuous(name ="Percentage of trees with Porifera-sister", breaks = seq(0,100,10), limits = c(0,110)) +
  scale_x_continuous(name = "Number of sites / Proportion of sites") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 20, margin = margin(t = 15, r = 0, b = 10, l = 0)),
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 15, b = 0, l = 10)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        strip.text = element_text(size = 20))
p_path <- paste0(plot_dir, "exploratory_NumSites_Porifera-sister.")
ggsave(filename = paste0(p_path, "png"), plot = p, device = "png")
ggsave(filename = paste0(p_path, "pdf"), plot = p, device = "pdf")

### Plot 3: Percent of trees with Cten/Cnid paraphyletic against number of sites ###
plot_df <- melt(al_df, 
                id.vars = c("ID", "best_model", "num_taxa", "num_sites", "percent_CTEN_CNID_monophyletic"),
                measure.vars = c("num_sites", "proportion_invariant_sites", "proportion_informative_sites") )
p <- ggplot(data = plot_df, aes(x = value, y = percent_CTEN_CNID_monophyletic)) +
  facet_grid(.~variable, scales = "free_x", labeller = labeller(variable = var_labs)) +
  geom_point(size = 3, color = "darkblue", alpha = 0.4) +
  scale_y_continuous(name ="Percentage of trees with\nmonophyletic (Ctenophora+Cnidaria)", breaks = seq(0,100,10), limits = c(0,110)) +
  scale_x_continuous(name = "Number of sites / Proportion of sites") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 20, margin = margin(t = 15, r = 0, b = 10, l = 0)),
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 15, b = 0, l = 10)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        strip.text = element_text(size = 20))
p_path <- paste0(plot_dir, "exploratory_NumSites_CtenophoraCnidaria-monophyletic.")
ggsave(filename = paste0(p_path, "png"), plot = p, device = "png")
ggsave(filename = paste0(p_path, "pdf"), plot = p, device = "pdf")

### Plot 4: Percent of trees with Ctenophora sister against length of branch to Ctenophora or Porifera clades ###
# Create labeller function
var_labs <- c("Ctenophora", "Porifera")
names(var_labs) = c("ctenophora_clade_branch_length", "porifera_clade_branch_length")
# Create palette
var_palette <- c(metazoan_clade_palette["Ctenophora"], metazoan_clade_palette["Porifera"])
names(var_palette) <- c("ctenophora_clade_branch_length", "porifera_clade_branch_length")
# Create plot
plot_df <- melt(al_df, 
                id.vars = c("ID", "best_model", "num_taxa", "num_sites", "percent_CTEN_sister", "percent_PORI_monophyletic", "percent_CTEN_CNID_monophyletic"),
                measure.vars = c("ctenophora_clade_branch_length", "porifera_clade_branch_length") )
p <- ggplot(data = plot_df, aes(x = value, y = percent_CTEN_sister, color = variable)) +
  facet_wrap(variable~., labeller = labeller(variable = var_labs)) +
  geom_point(size = 3, alpha = 0.6) +
  labs(title = "Length of branch leading to Metazoan clades") +
  scale_y_continuous(name ="Percentage of trees with Ctenopora-sister", breaks = seq(0,100,10), limits = c(0,100)) +
  scale_x_continuous(name = "Branch length (subs/site)") +
  scale_colour_manual(values = var_palette, guide = "none") +
  theme_bw() +
  theme(plot.title = element_text(size = 25, hjust = 0.5, margin = margin(t = 15, r = 0, b = 20, l = 0)),
        axis.title.x = element_text(size = 20, margin = margin(t = 15, r = 0, b = 10, l = 0)),
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 15, b = 0, l = 10)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        strip.text = element_text(size = 20))
p_path <- paste0(plot_dir, "exploratory_BranchLength_Ctenophora-sister.")
ggsave(filename = paste0(p_path, "png"), plot = p, device = "png")
ggsave(filename = paste0(p_path, "pdf"), plot = p, device = "pdf")

### Plot 5: Percent of trees with Ctenophora sister against clade depth for Ctenophora or Porifera clades ###
# Create labeller function
var_labs <- c("Ctenophora", "Porifera")
names(var_labs) = c("ctenophora_clade_depth", "porifera_clade_depth")
# Create palette
var_palette <- c(metazoan_clade_palette["Ctenophora"], metazoan_clade_palette["Porifera"])
names(var_palette) <- c("ctenophora_clade_depth", "porifera_clade_depth")
# Create plot
plot_df <- melt(al_df, 
                id.vars = c("ID", "best_model", "num_taxa", "num_sites", "percent_CTEN_sister", "percent_PORI_monophyletic", "percent_CTEN_CNID_monophyletic"),
                measure.vars = c("ctenophora_clade_depth", "porifera_clade_depth") )
p <- ggplot(data = plot_df, aes(x = value, y = percent_CTEN_sister, color = variable)) +
  facet_wrap(variable~., labeller = labeller(variable = var_labs)) +
  geom_point(size = 3, alpha = 0.6) +
  labs(title = "Depth of Metazoan clades") +
  scale_y_continuous(name ="Percentage of trees with Ctenopora-sister", breaks = seq(0,100,10), limits = c(0,100)) +
  scale_x_continuous(name = "Clade depth (subs/site)\nCalculated from max. branching.times") +
  scale_colour_manual(values = var_palette, guide = "none") +
  theme_bw() +
  theme(plot.title = element_text(size = 25, hjust = 0.5, margin = margin(t = 15, r = 0, b = 20, l = 0)),
        axis.title.x = element_text(size = 20, margin = margin(t = 15, r = 0, b = 10, l = 0)),
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 15, b = 0, l = 10)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        strip.text = element_text(size = 20))
p_path <- paste0(plot_dir, "exploratory_CladeDepth_Ctenophora-sister.")
ggsave(filename = paste0(p_path, "png"), plot = p, device = "png")
ggsave(filename = paste0(p_path, "pdf"), plot = p, device = "pdf")

