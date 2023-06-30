## caitlinch/metazoan-mixtures/code/func_plotting.R
# Functions for plotting graphs and trees
# Caitlin Cherryh 2023


#### Packages ####
library(ggtree)


#### Functions to plot figures for introduction ####
bw.monophyletic.clades.plot <- function(tree){
  # Function to plot in black and white overview of hypotheses for relationships between 
  #   clades of the Animal tree of life, where the sponge clade is monophyletic
  
  # Create the plot
  temp_plot <- ggtree(tree) + geom_tree() + 
    geom_tiplab(geom = "text", size = 10) +
    theme_tree2(plot.margin=margin(6, 145, 6, 6)) +
    coord_cartesian(clip = "off") +
    theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"))
  # Return the plot
  return(temp_plot)
}

bw.paraphyletic.clades.plot <- function(tree, label_nodes = c(NA,NA)){
  # Function to plot in black and white overview of hypotheses for relationships between 
  #   clades of the Animal tree of life, where the sponge clade is paraphyletic
  
  # Create the plot
  temp_plot <- ggtree(tree) + geom_tree() + 
    geom_tiplab(geom = "text", size = 8) +
    theme_tree2(plot.margin=margin(6, 70, 6, 6)) +
    coord_cartesian(clip = "off") +
    geom_cladelab(node=label_nodes[1], label="Porifera", textcolor="grey50", barcolor = "grey50", offset=7, geom = "text", align=TRUE, fontsize = 7, hjust = -0.1) +
    geom_cladelab(node=label_nodes[2], label="Porifera", textcolor="grey50", barcolor = "grey50", offset=7, geom = "text", align=TRUE, fontsize = 7, hjust = -0.1) +
    theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"))
  # Return the plot
  return(temp_plot)
}

color.clades.plot <- function(tree, tip_labels, color_palette, save.plot = FALSE, output_directory, output_id){
  # Function to plot in black and white overview of hypotheses for relationships between 
  #   clades of the Animal tree of life, where the sponge clade is monophyletic
  
  # Create the plot
  temp_plot <- ggtree(tree) %<+% tip_labels +
    geom_tree(size = 1.5) +
    geom_tiplab(aes(label = lab, color = color), parse = TRUE, show.legend = FALSE, geom = "text", size = 10, fontface = 2) +
    theme_tree2(plot.margin=margin(6, 145, 6, 6)) +
    coord_cartesian(clip = "off") +
    theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white")) +
    scale_color_manual(values = color_palette)
  # Save the plot if desired
  if (save.plot == TRUE){
    # Set output parameters
    if (grepl("1|2|3", "hypothesis_tree_1_plot_color")){
      # For hypothesis trees 1,2,3 (monophyletic porifera)
      plot_params <- c("png_width" = 974, "png_height" = 723, "svg_width" = 5, "svg_height" = 3.5)
    } else if (grepl("4|5", "hypothesis_tree_1_plot_color")){
      # For hypothesis trees 4,5 (paraphyletic porifera)
      plot_params <- c("png_width" = 1364, "png_height" = 723, "svg_width" = 7, "svg_height" = 3.5)
    }
    # Output png
    hypothesis_plot_file <- paste0(output_directory, output_id, ".png")
    png(filename = hypothesis_plot_file, width = plot_params[["png_width"]], height = plot_params[["png_height"]], units = "px", pointsize = 12, bg = "white")
    temp_plot
    dev.off()
    # Output svg
    hypothesis_plot_file <- paste0(output_directory, output_id, ".svg")
    svg(filename = hypothesis_plot_file, width = plot_params[["svg_width"]], height = plot_params[["svg_height"]], bg = "white")
    temp_plot
    dev.off()
  }
  # Return the plot
  return(temp_plot)
}




