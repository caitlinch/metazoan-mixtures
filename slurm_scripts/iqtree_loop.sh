#!/bin/bash

# Change to folder where you want to create the files
cd /Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/ML_iqtree_sh_scripts/

# Split the iqtree calls into 18 files of 24 lines each
split -l 24 --numeric-suffixes[=1] /Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/maximum_likelihood_iqtree2_calls.csv split_

# Add file extensions to split_ files
# Add .txt to all filenames...unless they are already .txt
for f in ls split_*; do case "$f" in *.txt) echo skipped $f;; *) mv "$f" "$f".txt; esac; done

# Create a .sh file for each tree
for i in {1..18} 
do
	# Copy the template file to a new file where the name includes the row number
	cp /Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/slurm_scripts/iqtree_format.sh iqtree_run_$i.sh
	# Update the row number in the new file for the sbatch commands
	sed -i'.bak' "s/\${id}/$i/" iqtree_run_$i.sh
	# Append the split_$i file to add the iqtree command lines for this file
	cat split_$i.txt >> iqtree_run_$i.sh
done

# Remove all files with .bak extension
ls *.bak
rm *.bak

# Create a list of all the .sh files
ls iqtree_run_* > /Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/iqtree_individual_job_runs.sh

# Add the prefix "sbatch " to each line (to call all  432 jobs at once)
sed -i'.bak' 's/^/sbatch /' /Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/iqtree_individual_job_runs.sh


