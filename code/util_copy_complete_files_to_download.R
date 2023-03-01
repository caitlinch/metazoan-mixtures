# Copy complete files to download

# Establish file paths
old_dir <- "/home/u5348329/metazoan-mixtures/output/maximum_likelihood_trees/"
new_dir <- "/home/u5348329/metazoan-mixtures/output/ml_tree_output_files/"

# Open parameters dataframe
df <- read.csv("/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/output/maximum_likelihood_tree_estimation_parameters.csv", stringsAsFactors = FALSE)
all_iqfiles <- df$iqtree_file
# Get the list of which iqtree files exists
existing_iqfiles <- all_iqfiles[which(file.exists(paste0(old_dir, all_iqfiles)) == TRUE)]

# Copy all files that do not exist in the new directory
for (iqf in existing_iqfiles){
  if (file.exists(paste0(new_dir, iqf)) == FALSE){
    file.copy(from = paste0(old_dir, iqf), to = paste0(new_dir, iqf), copy.date = T)
  }
  if (file.exists(paste0(new_dir, gsub("\\.iqtree", "\\.treefile", iqf))) == FALSE){
    file.copy(from = paste0(old_dir, gsub("\\.iqtree", "\\.treefile", iqf)), to = paste0(new_dir, gsub("\\.iqtree", "\\.treefile", iqf)), copy.date = T)
  }
  if (file.exists(paste0(new_dir, gsub("\\.iqtree", "\\.log", iqf))) == FALSE){
    file.copy(from = paste0(old_dir, gsub("\\.iqtree", "\\.log", iqf)), to = paste0(new_dir, gsub("\\.iqtree", "\\.log", iqf)), copy.date = T)
  }
}
