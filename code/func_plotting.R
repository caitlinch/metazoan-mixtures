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

color.clades.plot <- function(tree, tip_labels, color_palette, xlimits){
  ## Function to plot in black and white overview of hypotheses for relationships between 
  ##   clades of the Animal tree of life, where the sponge clade is monophyletic
  # Check whether x limits are provided
  check_xlims <- is.numeric(xlimits)
  # Create plot
  if (check_xlims == FALSE){
    # Create the plot without specifying x limits)
    temp_plot <- ggtree(tree) %<+% tip_labels +
      geom_tree(size = 1.5) +
      geom_tiplab(aes(label = lab, color = color), parse = TRUE, show.legend = FALSE, geom = "text", size = 10, fontface = 2) +
      theme_tree2(plot.margin=margin(6, 145, 6, 6)) +
      coord_cartesian(clip = "off") +
      theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white")) +
      scale_color_manual(values = color_palette)
  } else if (check_xlims == TRUE){
    # Create the plot with specified x limits
    temp_plot <- ggtree(tree) %<+% tip_labels +
      geom_tree(size = 1.5) +
      geom_tiplab(aes(label = lab, color = color), parse = TRUE, show.legend = FALSE, geom = "text", size = 10, fontface = 2) +
      xlim(xlimits[[1]], xlimits[[2]]) +
      theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white")) +
      scale_color_manual(values = color_palette)
  }
    # Return the plot
    return(temp_plot)
}




