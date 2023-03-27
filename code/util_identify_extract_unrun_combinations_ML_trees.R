# Check which trees still need running

# Set directories
ml_dir <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/output/maximum_likelihood_trees/"

# Open parameters dataframe
df <- read.table("/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/output/01_01_maximum_likelihood_tree_estimation_parameters.tsv", header = TRUE, stringsAsFactors = FALSE)
all_iqfiles <- df$iqtree_file
# Get all treefiles from the output folder
all_files <- list.files(ml_dir)
iqfiles <- grep("\\.iqtree", all_files, value = T)
# Get the list of completed files
completed_files <- iqfiles[file.exists(paste0(ml_dir, iqfiles))]
incomplete_files <-  df$iqtree_file[setdiff(1:nrow(df), which(completed_files %in% df$iqtree_file))]
incomplete_prefix = df$prefix[setdiff(1:nrow(df), which(completed_files %in% df$iqtree_file))]
# Create dataframe of files to run
df2 <- data.frame(index = 1:length(incomplete_files), prefix = incomplete_prefix, file = incomplete_files)
# Save output dataframe
df2_file_path <- paste0("/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/output/missing_runs_", gsub("-","", Sys.Date()),".csv")
write.csv(df2, file = df2_file_path, row.names = F)
