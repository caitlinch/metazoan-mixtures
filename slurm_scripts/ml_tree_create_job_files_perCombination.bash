#!/bin/bash

# Change to folder where you want to create the files
# cd /Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/ML_iqtree_sh_scripts/
cd /home/u5348329/metazoan-mixtures/slurm_files/ML_tree_iqtree_runs_individual/

# Create a .sh file for each tree
for i in {1..384} 
do
	# Copy the template file to a new file where the name includes the row number
	# cp /Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/slurm_scripts/iqtree_format.sh iqtree_run_$i.sh
	cp ../job_file_template.sh iqtree_run_$i.sh
	# Update the row number in the new file for the sbatch commands
	sed -i'.bak' "s/\${id}/$i/" iqtree_run_$i.sh
	# Get the i^th line of the iqtree calls text file and append it to the iqtree run file for this id
    sed "${i}q;d" /home/u5348329/metazoan-mixtures/output/maximum_likelihood_iqtree2_calls.txt >> iqtree_run_$i.sh

done

# Remove all files with .bak extension
ls *.bak
rm *.bak

# Create a list of all the .sh files
ls iqtree_run_* > iqtree_individual_job_runs.sh

# Add the prefix "sbatch " to each line (to call all 384 jobs at once)
sed -i'.bak' 's/^/sbatch /' iqtree_individual_job_runs.sh
# Remove backup file 
rm iqtree_individual_job_runs.sh.bak

# Convert file line endings from dos format to unix format using dos2unix
for i in {1..384}
do 
	dos2unix iqtree_run_$i.sh
done