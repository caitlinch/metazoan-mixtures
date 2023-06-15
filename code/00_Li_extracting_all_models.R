## caitlinch/metazoan-mixtures/code/00_matrix_parameters.R
# This script extracts the list of models used in a previous study, from the github repository of that study
# Caitlin Cherry, 2023


#### 1. Input parameters ####
## Specify parameters:
# supp_dir <- directory containing supplementary data from Li et. al. (2021) paper (downloaded from github repository https://github.com/dunnlab/animal_tree_root)
# main_dir <- directory for the GitHub repository for this project: caitlinch/metazoan-mixtures

supp_dir <- "/Users/caitlincherryh/Documents/Repositories_notMine/animal_tree_root-main" 
main_dir <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"



#### 2. Prepare variables, packages, and functions
# Source functions for extracting models and processing datasets
source(paste0(main_dir, "code/func_data_processing.R"))



#### 3. Extract vector of all models used in Li et. al. (2021) ####
# List all files in the supp_dir
all_files <- list.files(supp_dir, recursive = TRUE, full.names = TRUE)
# Extract files that contain information about models
files_of_interest <- c(grep("Supplementary_Table_1.csv", all_files, value = TRUE),
                       grep("Supplementary_Table_2.csv", all_files, value = TRUE))

### Extract the list of models from each relevant supplementary table:
## Supplementary Table 1
#     Summary of a total of 164 phylogenomic analyses were transcribed from the literature
#     (Table is converted from analyses_published in Rdata).
st1 <- read.csv(files_of_interest[[1]])
# Reduce to only AA models (no recoding)
st1 <- st1[st1$recoding == "aa", ]
# Remove partitioning models
st1 <- st1[st1$model_rate_matrix != "data partitioning", ]
# Remove NAs in the model_rate_matrix
st1$model_rate_matrix[is.na(st1$model_rate_matrix)] <- ""
# Relabel the model site rate heterogeneity column
st1$model_site_rate_hetero[is.na(st1$model_site_rate_hetero)] <- ""
st1$model_site_rate_hetero[st1$model_site_rate_hetero == "GAMMA"] <- "G"
# Relabel the model equilibrium column
st1$model_equilibrium[st1$model_equilibrium == "mix"] <- ""
st1$model_equilibrium[is.na(st1$model_equilibrium)] <- ""
st1$model_equilibrium[st1$model_equilibrium == "Empritical"] <- "F"
st1$model_equilibrium[st1$model_equilibrium == "Estimated"] <- "FO"
# Split into two sections
st1_f <- st1[st1$model_equilibrium == "F" | st1$model_equilibrium == "FO",]
st1_no.f <- st1[st1$model_equilibrium != "F" & st1$model_equilibrium != "FO",]
# Paste models together
st1_models_1 <- paste0(st1_f$model_rate_matrix, "+", st1_f$model_site_rate_hetero, "+", st1_f$model_equilibrium)
st1_models_2  <- paste0(st1_no.f$model_rate_matrix, "+", st1_no.f$model_equilibrium, "+", st1_no.f$model_site_rate_hetero)
# Combine models
st1_models <- c(st1_models_1, st1_models_2)
# Remove any pluses from the models
st1_models <- unlist(lapply(st1_models, remove.extra.plusses))
# Reduce to only the unique models
st1_models <- sort(unique(st1_models))

## Supplementary Table 2
#     Summary of a total of 106 phylogenomic analyses conducted in this study 
#     (Table is converted from analyses_new in R data).
st2 <- read.csv(files_of_interest[[2]])
# Extract list of models 
st2_models <- st2$model
# Reduce to only the unique models
st2_models <- sort(unique(st2_models))

# Combine and collate into single vector
li_models <- sort(unique(c(st1_models, st2_models)))
# Alphabetise models by model chunks
li_models <- unique(sort(unlist(lapply(li_models, sort.model.chunks))))

### Save the output
# Ensure the output folder exists (and create it if it doesn't)
output_dir <- paste0(main_dir, "data/")
if (dir.exists(output_dir) == FALSE){
  dir.create(output_dir)
}

# Save the models as a text file
output_path <- paste0(output_dir, "Li2021_all_models.txt")
write(li_models, file = output_path)


