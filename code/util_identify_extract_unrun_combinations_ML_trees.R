# Check which trees still need running

# Open parameters dataframe
df <- read.csv("/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/output/maximum_likelihood_tree_estimation_parameters.csv", stringsAsFactors = FALSE)
all_iqfiles <- df$iqtree_file
# Get all treefiles from the output folder
all_files <- list.files("/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/output/maximum_likelihood_trees/")
iqfiles <- grep("\\.iqtree", all_files, value = T)
# Create dataframe of files to run
inds_to_run <- setdiff(1:384, which(all_iqfiles %in% iqfiles))
ids_to_run <- df$prefix[inds_to_run]
files_to_run <- paste0("iqtree_run_", inds_to_run, ".sh")
df2 <- data.frame(index = inds_to_run, prefix = ids_to_run, sbatch_file = files_to_run)
# Save output dataframe
df2_file_path <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/output/missing_runs_20230309.csv"
write.csv(df2, file = df2_file_path, row.names = F)
