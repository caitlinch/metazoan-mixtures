# Input parameters
# supp_dir <- directory containing supplementary data from Li et. al. (2021) paper (downloaded from github repository)

supp_dir <- "/Users/caitlin/Downloads/animal_tree_root-main" 
mandir <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/animal_tree_root-main/"



#### 1. Prepare variables, packages, and functions
source(paste0(main_dir, "code/func_data_processing.R"))



#### 2. Extract vector of all models used in Li et. al. (2021) ####
# List all files in the supp_dir
all_files <- list.files(supp_dir, recursive = TRUE, full.names = TRUE)
# Extract files that contain information about models
files_of_interest <- c(grep("Supplementary_Table_1.csv", all_files, value = TRUE),
                       grep("Supplementary_Table_2.csv", all_files, value = TRUE),
                       grep("Supplementary_Table_5.csv", all_files, value = TRUE))

## Look at each table in turn
# Supplementary Table 1
st1 <- read.csv(files_of_interest[[2]])
# Reduce to only AA models (no recoding)
st1 <- st1[st1$recoding == "aa", ]
# Remove partitioning models
st1 <- st1[st1$model_rate_matrix != "data partitioning", ]
# Relabel the model site rate heterogeneity column
st1$model_site_rate_hetero[is.na(st1$model_site_rate_hetero)] <- ""
st1$model_site_rate_hetero[st1$model_site_rate_hetero == "GAMMA"] <- "G"
# Relabel the model equilibrium column
st1 <- st1[st1$model_equilibrium != "mix", ]
st1$model_equilibrium[is.na(st1$model_equilibrium)] <- ""
st1$model_equilibrium[st1$model_equilibrium == "Empritical"] <- "F"
# Split into two sections
st1_f <- st1[st1$model_equilibrium == "F",]
st1_no.f <- st1[st1$model_equilibrium != "F",]
# Paste models together
st1_models_1 <- paste0(st1_f$model_rate_matrix, "+", st1_f$model_site_rate_hetero, "+", st1_f$model_equilibrium)
st1_models_2  <- paste0(st1_no.f$model_equilibrium, "+", st1_no.f$model_rate_matrix, "+", st1_no.f$model_site_rate_hetero)
# Combine models
st1_models <- c(st1_models_1, st1_models_2)
# Remove any pluses from the models
st1_models <- unlist(lapply(st1_models, remove.extra.pluses))
# Reduce to only the unique models
st1_models <- sort(unique(st1_models))

# Supplementary Table 2
st2 <- read.csv(files_of_interest[[3]])
# Extract list of models 
st2_models <- st2$model
# Reduce to only the unique models
st2_models <- unique(st2_models)

# Supplementary Table 5
st5 <- read.csv(files_of_interest[[4]])
