## caitlinch/metazoan-mixtures/code/func_data_processing.R
# Functions for manipulating, processing, and preparing datasets
# Caitlin Cherryh 2023


#### Packages ####
library(stringr) # Used in extract.model.details for detecting presence/absence of a pattern in a string (using str_detect) 
library(ape) # Used in check.tree.taxa to read in a phylogenetic tree from file (using read.tree)



#### Process models to prepare for IQ-Tree runs ####
create.model.dataframe <- function(model_vector, output.counts = FALSE, output.model.chunks = FALSE){
  # Function to read in a partition file and return the list of models applied to the charsets
  
  # Identify how many times each model was used
  model_t <- table(model_vector)
  model_count <- as.numeric(model_t)
  model_name <- names(model_t)
  
  # Construct a dataframe with the model name (and counts if required)
  if (output.counts == TRUE){
    model_df <- data.frame(model = model_name,
                           count = model_count)
  } else if (output.counts == FALSE){
    model_df <- data.frame(model = model_name)
  }
  # Get the starting names from the model_df
  model_df_init_names <- names(model_df)
  
  # If desired, output columns for the model chunks
  if (output.model.chunks == TRUE){
    # Find the maximum number of chunks that belong to each model
    max_number_chunks <- max(unlist(lapply(strsplit(model_name, "\\+"), length)))
    # Break each model down into the component parts
    for (i in 1:max_number_chunks){
      # Extract the chunks at this position
      model_chunks <- unlist(lapply(model_name, get.model.chunk, i))
      # Add the chunks at this position to the model dataframe
      model_df <- cbind(model_df, model_chunks)
      # Create a column name for this position
      col_names <- c(model_df_init_names, paste0("model_position_",1:i))
      # Rename the columns for the model dataframe
      names(model_df) <- col_names
    } # end for (i in 1:max_number_chunks)
  } # end if (output.model.chunks == TRUE)
  
  # Check whether each model contains +I, +F, or +G. Add results to model dataframe
  model_df$plus_f <- grepl("\\+F", model_name)
  model_df$plus_i <- grepl("\\+I\\+", model_name)
  model_df$plus_g <- grepl("\\+G", model_name)
  
  # Return model dataframe
  return(model_df)
}


extract.partition.models <- function(partition_file){
  # Function to read in a partition file and return the list of models applied to the charsets
  
  # Open the partition file
  p <- readLines(partition_file)
  # Find the "charpartition" section (this lists the model for each defined charset)
  start_ind <- grep("charpartition", p, ignore.case = T)
  # Find the end of the charpartition section (find the line containing "end" that is after the start of the charpartition section)
  end_ind <- grep("end", p, ignore.case = T)[which(grep("end", p, ignore.case = T) > start_ind)]
  # Extract all lines in the charpartition section
  charpartition_lines <- p[start_ind:end_ind]
  
  # Identify all lines containing a colon (":") - these will specify a model for a charset
  model_lines <- grep(":", charpartition_lines, value = T)
  # Split the lines at the colon and select the first section of each line
  models <- unlist(lapply(model_lines, 
                          function(l){gsub(" ", "", strsplit(l, ":")[[1]][[1]])} ))
  
  # Return the list of models (one model per charset)
  return(models)
}


get.model.chunk <- function(model, number){
  # Small function to break up a model into chunks and return the chunk at number position (e.g. first position, third position)
  
  # Break model into chunks
  m_chunks <- strsplit(model, "\\+")[[1]]
  # Identify how many chunks there are
  number_chunks <- length(m_chunks)
  
  # Select the chunk at position number
  if (number <= number_chunks){
    # If the selected position is present within the model, extract it
    chunk <- m_chunks[number]
  } else if (number > number_chunks){
    # If the number is higher than the number of chunks within the model, return NA
    chunk <- NA
  }
  
  # Return the chunk at position number for this model
  return(chunk)
}


remove.extra.plusses <- function(m){
  # Quick function to remove any extra pluses in a model of sequence evolution
  
  # Change any double plusses ("++") into a single plus ("+)
  m <- gsub("\\+\\+", "+", m)
  
  # Remove any pluses at the start or the end of the model
  # Split the model at the plus sign
  m_split <- strsplit(m, "\\+")[[1]]
  # Remove any empty objects (i.e. "")
  m_split <- m_split[m_split != ""]
  # If there are any chunks in the model (i.e. if the model was not a + or ++,),
  # join together the chunks from the model
  if (length(m_split) > 0){
    # Paste the model back together using a plus sign
    new_m <- paste(m_split, collapse = "+")
    
    if (new_m != "" & new_m != "+"){
      # Do not return new_m if it is empty ("") or a single plus ("+")
      return(new_m)
    } # end if (new_m != "" | new_m != "+")
  } # end if (length(m_split) > 0)
} # end remove.extra.plusses <- function(m)


sort.model.chunks <- function(m){
  # Small function to sort the chunks in a model and return it in a standard order
  
  # Break model into chunks
  m_chunks <- strsplit(m, "\\+")[[1]]
  # Remove any special chunks
  # Special chunks: +I, +G, +G4, +FO, +F, +R, +R4
  main_chunks <- m_chunks[which(!(m_chunks %in% c("FO", "F", "G", "G4", "I", "R", "R4")))]
  # Identify the special chunks in the model
  special_chunks <- m_chunks[which((m_chunks %in% c("FO", "F", "G", "G4", "I", "R", "R4")))]
  # Sort the list of model chunks
  new_m_vec <- c(sort(main_chunks), special_chunks)
  # Stick the model chunks together
  new_m <- paste(new_m_vec, collapse = "+")
  # Return the new, sorted model
  return(new_m)
}




#### Extract details from IQ-Tree output files ####
extract.best.model <- function(iqtree_file){
  # Function that will extract the best model of sequence evolution or the model of sequence evolution used,
  #   given a .iqtree file
  
  # Check if the file is a PMSF file
  # The PMSF models do not include a model search or a BIC score for the model
  pmsf_check <- grepl("PMSF", iqtree_file)
  if (pmsf_check == FALSE){
    # The model is not a PMSF model. Continue to extract best model
    # Check if the iqtree file exists
    if (file.exists(iqtree_file) == TRUE){
      # If the iqtree_file does exist:
      ## Open the .iqtree file:
      iq_lines <- readLines(iqtree_file)
      
      ## Check for a ModelFinder section:
      # Determine whether there is a ModelFinder section
      mf_ind <- grep("ModelFinder", iq_lines)
      # Determine whether there is a line detailing the best model
      bm_ind <- grep("Best-fit model according to", iq_lines)
      
      ## Check for a Substitution Process section:
      # Determine the starting line of this section
      sp_ind <- grep("SUBSTITUTION PROCESS", iq_lines)
      # Determine the line detailing the model used
      mos_ind <- grep("Model of substitution", iq_lines)
      
      ## Extract the best fit model from the .iqtree file:
      if ((identical(mf_ind, integer(0)) == FALSE) & (identical(bm_ind, integer(0)) == FALSE)){
        # If ModelFinder was run, extract the best model from the ModelFinder section of the .iqtree file
        # Extract the line containing the best fit model
        m_line <- iq_lines[bm_ind]
      } else if ((identical(sp_ind, integer(0)) == FALSE) & (identical(mos_ind, integer(0)) == FALSE)) {
        # If there is no ModelFinder section, extract the model used from the substitution process section
        m_line <- iq_lines[mos_ind]
      } else {
        m_line <- "NA:NA"
      }
      
      ## Format the model nicely for output: 
      # Split the line at the colon into two parts
      m_line_split <- strsplit(m_line, ":")[[1]]
      # If the best model is a single model, the length of m_line_split will be 2
      #     One section for the explanatory text and one for the model
      # If the best model is a partition model, it will have more than two sections when split by colons
      # Extract the second part of the line onwards (contains the best fit model)
      best_model <- m_line_split[2:length(m_line_split)]
      # If best_model is longer than 1, paste it together again using colons
      if (length(best_model) >1){
        best_model <- paste(best_model, collapse = ":")
      }
      # Remove any white space from the best model
      best_model <- gsub(" ", "", best_model)
    } else if (file.exists(iqtree_file) == FALSE){
      # If the iqtree_file doesn't exist, return NA
      best_model = NA
    } # end if (file.exists(iqtree_file) == TRUE){
  } else if (pmsf_check == TRUE){
    # If the PMSF model was used, there is no modelfinder output to search and no comparison BIC
    # Return NA
    best_model = NA
  } # end if (pmsf_check == FALSE){
  # Return the best model from the iqtree_file (if the file exists)
  return(best_model)
}


extract.treefile <- function(tree_file){
  # Small function to extract a tree using a file path
  
  # Check whether the tree file exists
  if (file.exists(tree_file) == TRUE){
    # If the file exists, open and read the tree
    tree_text <- readLines(tree_file)
  } else if (file.exists(tree_file) == FALSE){
    # If the file doesn't exist, return NA
    tree_text <- NA
  }
  
  # Return the tree (if it exists)
  return(tree_text)
}  


extract.branch.length.wrapper <- function(row_id, alignment_df, tree_directory, clade){
  ## Wrapper function for extract.branch.length
  
  # Extract row information
  row <- alignment_df[row_id,]
  # Call function
  bl <- extract.branch.length(dataset = row$dataset, matrix_name = row$matrix_name, best_model = row$best_model,
                              tree_directory = tree_directory, clade = clade)
  # Return branch length
  return(bl)
}


extract.branch.length <- function(dataset, matrix_name, best_model, tree_directory, clade){
  ## Function to extract branch length for either Ctenophora or Porifera clade
  # Find and open best tree for this dataset
  all_trees <- list.files(tree_directory)
  row_tree_file <- grep(best_model, grep(matrix_name, grep(dataset, all_trees, value = T), value = T), value = T)
  row_tree_file_path <- paste0(tree_directory, row_tree_file)
  raw_tree <- read.tree(row_tree_file_path)
  # Extract clades from tip labels
  outgroup_species <- grep("outgroup", raw_tree$tip.label, value = T, ignore.case = T)
  ctenophora_species <- grep("ctenophora", raw_tree$tip.label, value = T, ignore.case = T)
  porifera_species <- grep("porifera", raw_tree$tip.label, value = T, ignore.case = T)
  # Root at outgroup
  og_tips <- raw_tree$tip.label
  tree <- root(raw_tree, outgroup = outgroup_species, resolve.root = T)
  # Get node, branch numbers and branch lengths for Ctenophora clade
  if (length(ctenophora_species) > 1){
    # If multiple sponge species
    ctenophora_node <- getMRCA(tree, ctenophora_species)
    ctenophora_branch <- which(tree$edge[,2] == ctenophora_node)
    ctenophora_branch_length <- tree$edge.length[ctenophora_branch]
  } else {
    # If single sponge species 
    ctenophora_branch <- which(tree$edge[,2] == which(tree$tip.label == ctenophora_species))
    ctenophora_node <- tree$edge[ctenophora_branch, 1]
    ctenophora_branch_length <- NA
  }
  # Get node, branch numbers and branch lengths for Porifera clade
  if (length(porifera_species) > 1){
    # If multiple sponge species
    porifera_node <- getMRCA(tree, porifera_species)
    porifera_branch <- which(tree$edge[,2] == porifera_node)
    porifera_branch_length <- tree$edge.length[porifera_branch]
  } else {
    # If single sponge species 
    porifera_branch <- which(tree$edge[,2] == which(tree$tip.label == porifera_species))
    porifera_node <- tree$edge[porifera_branch, 1]
    porifera_branch_length <- NA
  }
  # Return requested branch lengths
  if (clade == "Porifera"){
    return(porifera_branch_length)
  } else if (clade == "Ctenophora"){
    return(ctenophora_branch_length)
  }
}


extract.model.log.likelihood <- function(iqtree_file, var = "LogL"){
  # Function to extract the log likelihood of a model from in IQ-Tree
  # Can extract either the log likelihood (var = "LogL), the BIC (var = "BIC),
  #   the weighted BIC (var = "wBIC") or all three (var = "All")
  
  
  # Check if the file is a PMSF file
  # The PMSF models do not include a model search or a BIC score for the model
  pmsf_check <- grepl("PMSF", iqtree_file)
  if (pmsf_check == FALSE){
    # The model is not a PMSF model. Continue to extract model parameters
    # Check if the iqtree file exists
    if (file.exists(iqtree_file) == TRUE){
      # If the iqtree_file does exist:
      ## Open the .iqtree file:
      iq_lines <- readLines(iqtree_file)
      
      ## Extract the table of models sorted by BIC scores
      # Find the line detailing the best model
      table_ind <- grep("List of models sorted by BIC scores: ", iq_lines)
      # Check whether the table is present
      if (identical(table_ind, integer(0)) == TRUE){
        # No table of models present
        # Therefore no model log likelihood or BIC scores to return
        # Return NA for all outputs
        if (var == "All"){
          output = c("BestModel" = NA, "LogL" = NA, "BIC" = NA, "w-BIC" = NA)
        } else {
          output = NA
        } # end if (var == "All"){
      } else if (identical(table_ind, integer(0)) == FALSE){
        # Table of models is present
        # Add two to the table ind to get the start of the table
        table_start<- table_ind + 2
        table_end <- grep("AIC, w-AIC   : Akaike information criterion scores and weights.", iq_lines) - 2
        # Extract the table lines
        table_lines <- iq_lines[table_start:table_end]
        # Split table into rows using " "
        split_table_lines <- strsplit(table_lines, " ")
        # Remove empty strings from table lines
        split_table_lines <- lapply(1:length(split_table_lines), function(i){split_table_lines[[i]][split_table_lines[[i]] != ""]})
        # Remove "+" from table lines
        split_table_lines <- lapply(1:length(split_table_lines), function(i){split_table_lines[[i]][split_table_lines[[i]] != "+"]})
        # Remove "-" from table lines
        split_table_lines <- lapply(1:length(split_table_lines), function(i){split_table_lines[[i]][split_table_lines[[i]] != "-"]})
        # Make table lines into data frame
        table_df <- as.data.frame(do.call(rbind, split_table_lines[2:length(split_table_lines)]))
        names(table_df) <- split_table_lines[[1]]
        
        ## Extract the log likelihood, BIC and wBIC for the best model
        # Prepare output for return
        if (var == "LogL"){
          output = table_df$LogL[1]
        } else if (var == "BIC"){
          output = table_df$BIC[1]
        } else if (var == "wBIC"){
          output = table_df$`w-BIC`[1]
        } else if (var == "All"){
          output = c("BestModel" = table_df$Model[1], "LogL" = table_df$LogL[1], "BIC" = table_df$BIC[1], "w-BIC" = table_df$`w-BIC`[1])
        } # end if (var == "LogL")
      } # end if (identical(table_ind, integer(0)) == TRUE){
    } # end if (file.exists(iqtree_file) == TRUE)
  } else if (pmsf_check == TRUE){
    # If the PMSF model was used, there is no modelfinder output to search and no output value
    # Return NA
    # Prepare output for return
    if (var == "All"){
      # If output is named vector, return NA for each value in the vector
      output = c("BestModel" = NA, "LogL" = NA, "BIC" = NA, "w-BIC" = NA)
    } else {
      # If output is a single value - return NA
      output = NA
    } # end if (var == "All"){
  }# end if if (pmsf_check == FALSE){
  
  # Return the output
  return(output)
}


extract.tree.log.likelihood <- function(iqtree_file, var = "LogL"){
  # Function to extract the log likelihood of a model from in IQ-Tree
  # Can extract either the log likelihood (var = "LogL), the unconstrained log likelihood (without tree) (var = "ULL"),
  #   the number of free parameters (var = "NFP), the BIC (var = "BIC"), the total tree length (var = "TrLen"), 
  #   the sum of internal branch lengths (var = "SIBL"),  or all of the above (var = "All")
  
  # Check if the iqtree file exists
  if (file.exists(iqtree_file) == TRUE){
    # If the iqtree_file does exist:
    ## Open the .iqtree file:
    iq_lines <- readLines(iqtree_file)
    
    ## Extract the variables from the iqtree file
    if (var == "LogL"){
      # Log Likelihood
      output <- gsub(" ", "", strsplit(strsplit(iq_lines[grep("Log-likelihood of the tree\\:", iq_lines)], "\\:")[[1]][2], "\\(")[[1]][1])
    } else if (var == "ULL"){
      # Unconstrained log likelihood (without tree)
      output <- gsub(" ", "", strsplit(iq_lines[grep("Unconstrained log-likelihood \\(without tree\\)\\:", iq_lines)], "\\:")[[1]][2])
    } else if (var == "NFP"){
      # Number of free parameters
      output <- gsub(" ", "", strsplit(iq_lines[grep("Number of free parameters \\(\\#branches \\+ \\#model parameters\\)\\:", iq_lines)], "\\:")[[1]][2])
    } else if (var == "BIC"){
      # BIC
      output <- gsub(" ", "", strsplit(iq_lines[grep("Bayesian information criterion \\(BIC\\) score\\:", iq_lines)], "\\:")[[1]][2])
    } else if (var == "TrLen"){
      # Tree length
      output <- gsub(" ", "", strsplit(iq_lines[grep("Total tree length \\(sum of branch lengths\\)\\:", iq_lines)], "\\:")[[1]][2])
    } else if (var == "SIBL"){
      # Sum of internal branch lengths
      output <- gsub(" ", "", strsplit(strsplit(iq_lines[grep("Sum of internal branch lengths\\:", iq_lines)], "\\:")[[1]][2], "\\(")[[1]][1])
    } else if (var == "All"){
      # Extract all values
      logl <- gsub(" ", "", strsplit(strsplit(iq_lines[grep("Log-likelihood of the tree\\:", iq_lines)], "\\:")[[1]][2], "\\(")[[1]][1])
      ull <- gsub(" ", "", strsplit(iq_lines[grep("Unconstrained log-likelihood \\(without tree\\)\\:", iq_lines)], "\\:")[[1]][2])
      nfp <- gsub(" ", "", strsplit(iq_lines[grep("Number of free parameters \\(\\#branches \\+ \\#model parameters\\)\\:", iq_lines)], "\\:")[[1]][2])
      bic <- gsub(" ", "", strsplit(iq_lines[grep("Bayesian information criterion \\(BIC\\) score\\:", iq_lines)], "\\:")[[1]][2])
      treelen <- gsub(" ", "", strsplit(iq_lines[grep("Total tree length \\(sum of branch lengths\\)\\:", iq_lines)], "\\:")[[1]][2])
      sibl <- gsub(" ", "", strsplit(strsplit(iq_lines[grep("Sum of internal branch lengths\\:", iq_lines)], "\\:")[[1]][2], "\\(")[[1]][1])
      # Assemble all variables into a vector (for outputting all informating at once)
      output <- c("LogL" = logl, "Unconstrained_LogL" = ull, "NumFreeParams" = nfp, 
                  "BIC" = bic, "TreeLength" = treelen, "SumInternalBranchLengths" = sibl)
    } # end if (var == "LogL")
  } # end if (file.exists(iqtree_file) == TRUE)
  
  # Return the output
  return(output)
}


extract.rates <- function(iqtree_file){
  # Function to extract the rate parameters of a model from in IQ-Tree
  
  # Check if the iqtree file exists
  if (file.exists(iqtree_file) == TRUE){
    # If the iqtree_file does exist:
    ## Open the .iqtree file:
    iq_lines <- readLines(iqtree_file)
    
    # Check for +R parameters 
    rate_ind <- grep("Site proportion and rates\\:", iq_lines)
    # Check if the rate_ind is present (if yes, that means there's a line containing the rates and weights)
    if (identical(rate_ind,integer(0)) == FALSE & class(rate_ind) == "integer"){
      ## Rates (+R)
      # The site proportion and weights values are present
      # Extract the rate weights and parameters from the iqtree file
      rates_line <- iq_lines[rate_ind]
      # Remove text from the beginning of the line
      rates_raw <- strsplit(rates_line, "\\:")[[1]][2]
      # Split up by the spaces
      split_rates <- strsplit(rates_raw, " ")[[1]]
      # Remove any entries that are empty characters (i.e. "")
      split_rates <- split_rates[split_rates != ""]
      # Split rates again by the commas (",")
      split2_rates <- unlist(strsplit(split_rates, ","))
      # Remove brackets from values
      rate_vals <- unlist(strsplit(unlist(strsplit(split2_rates, "\\(")), "\\)"))
      # Make output (for attaching into IQ-Tree2 command line)
      # Order is: proportion_1, rate_1, proportion_2, rate_2, ....., proportion_n, rate_n
      rates_op <- paste(rate_vals, collapse = ",")
    } else {
      # The site proportion and weights values are missing
      # Return NA for this file
      rates_op <- NA
    } # end if (identical(rate_ind,logical(0)) == FALSE & class(rate_ind) == "integer")
  } # end if (file.exists(iqtree_file) == TRUE)
  
  # Return the output
  return(rates_op)
}


extract.gamma.values <- function(iqtree_file, gamma.parameter = "List"){
  # Function to extract the gamma parameters of a model from in IQ-Tree
  #   gamma.parameter = "List" - if model has gamma rate, returns the site relative weights and proportions
  #   gamma.parameter = "Shape" - if model has gamma rate, returns the gamma shape alpha 
  
  # Check if the iqtree file exists
  if (file.exists(iqtree_file) == TRUE){
    # If the iqtree_file does exist:
    ## Open the .iqtree file:
    iq_lines <- readLines(iqtree_file)
    
    # Check for +G parameters 
    gamma_shape_ind <- grep("Gamma shape alpha\\:", iq_lines)
    gamma_list_ind <- grep("Model of rate heterogeneity\\: Gamma", iq_lines)
    # Check if the gamma parameters are present (if yes, that means there's a line containing the rates and weights)
    if ( (identical(gamma_shape_ind,integer(0)) == FALSE & class(gamma_shape_ind) == "integer") |
         (identical(gamma_list_ind,integer(0)) == FALSE & class(gamma_list_ind) == "integer") ){
      ## Gamma (+G)
      if (gamma.parameter == "Shape"){
        ## Extract the value of the gamma shape alpha
        # Extract the line containing the gamma alpha value
        gamma_line <- iq_lines[gamma_shape_ind]
        # Process the line to extract just the gamma alpha value
        raw_gamma_alpha <- strsplit(gamma_line, "\\:")[[1]][2]
        gamma_alpha <- as.numeric(gsub(" ", "", raw_gamma_alpha))
        # Set the output 
        gamma_op <- gamma_alpha
      } else if (gamma.parameter == "List"){
        ## Extract the gamma weights and proportions
        #Create the matrix for discrete gamma categories
        g_start <- grep(" Category", iq_lines) + 1 # get the index for the first line of the gamma categories matrix
        empty   <- which(iq_lines == "") # get indexes of all empty lines
        empty   <- empty[empty > g_start] # get empty lines above gamma categories matrix
        g_end   <- empty[1] - 1 # get end index for gamma categories matrix (one less than next empty line)
        end_line <- iq_lines[g_end]
        # if the end isn't an empty line, subtract one from the end count 
        # to exclude lines like "Relative rates are computed as MEAN of the portion of the Gamma distribution falling in the category."
        # to see if this is what's happening, check whether the line starts with a numeric section (i.e. a category for the gamma rate)
        check_line <- length(strsplit(strsplit(end_line, "        " )[[1]][1], " ")[[1]])
        if (check_line > 3){
          # If the check_line is longer than 3 characters, it won't be a group for the gamma categories but an instruction
          # Instructions can be excluded from the gamma matrix (but categories can't)
          g_end = g_end - 1
        }
        # Extract the lines of interest
        g_lines <- iq_lines[g_start:g_end]
        # Split the g_lines at the spaces
        g_lines_split <- strsplit(g_lines, " ")
        # Remove empty values
        g_lines_neat <- lapply(1:length(g_lines_split), function(i){g_lines_split[[i]] <- g_lines_split[[i]][which(g_lines_split[[i]] != "")]})
        # Collect the values in the right order for the output
        raw_gamma_vals <- unlist(lapply(1:length(g_lines_neat), function(i){c(g_lines_neat[[i]][3], g_lines_neat[[i]][2])}))
        # Create a nice output format
        gamma_vals <- paste(gsub(" ", "", raw_gamma_vals), collapse = ",")
        # Set the output
        gamma_op <- gamma_vals
      } # end if (gamma.parameter == "Shape"){
    }else {
      # The site proportion and weights values are missing
      # Return NA for this file
      gamma_op <- NA
    } # end if ( (identical(gamma_shape_ind,integer(0)) == FALSE & class(gamma_shape_ind) == "integer") | ( ...) ){
  } # end if (file.exists(iqtree_file) == TRUE)
  
  # Return the output
  return(gamma_op)
}


extract.state.frequencies <- function(iqtree_file){
  # Given an iqtree file, this function will extract the state frequencies for the alignment
  
  ## Check if the iqtree file exists
  if (file.exists(iqtree_file) == TRUE){
    ## Open the .iqtree file:
    iq_lines <- readLines(iqtree_file)
    
    ## Extract state frequency details:
    # Check for presence of state frequencies details
    ind <- grep("State frequencies:", iq_lines)
    if (identical(ind, integer(0)) == FALSE){
      ## Determine whether state frequencies were empirically determined (or not!)
      # If there is a section for state frequencies, extract and output details 
      sf1 <- strsplit(iq_lines[[ind]], ":")[[1]][2]
      # Check whether state frequencies are needed
      sf1_squashed <- gsub(" ", "", sf1)
      if (sf1_squashed == "(empiricalcountsfromalignment)"){
        ## If state frequencies were determined from the alignment, 
        #    extract the lines with state frequencies from the .iqtree file
        # Get starting line for frequencies
        start_ind <- grep("State frequencies:", iq_lines) + 2
        # Take the 20 lines containing AA frequencies
        freq_lines <- iq_lines[start_ind:(start_ind+19)]
        # Split up the frequency lines into the label and the frequency
        freq_split <- unlist(strsplit(freq_lines, "="))
        
        ## Process the frequencies 
        # Get the frequency
        freq_nums <- freq_split[c(FALSE, TRUE)]
        # Remove any spaces (from IQTree formatting)
        freq_nums <- gsub(" ","",freq_nums)
        
        ## Process the labels
        # Get corresponding AA letter
        freq_names <- freq_split[c(TRUE, FALSE)]
        # Remove IQTree formatting
        freq_names <- gsub("pi\\(", "", freq_names)
        freq_names <- gsub("\\)", "", freq_names)
        freq_names <- gsub(" ", "", freq_names)
        
        ## Make the output pretty
        # Create a nice output by pasting together the frequencies in order
        f_op = paste(freq_nums, collapse = ",")
        
      } else if (sf1_squashed == "(equalfrequencies)"){
        ## If state frequencies were equal, generate equal state frequencies
        f_op <- paste(as.character(rep(1/20, 20)), collapse = ",")
      } else {
        f_op <- "State frequencies from model"
      }
    } else if (identical(ind, integer(0)) == TRUE){
      # If no details on state frequencies are needed, return an empty dataframe
      f_op <- NA
    }
  }
  
  ## Output the frequencies
  return(f_op)
}



extract.cat.frequencies <- function(iqtree_file){
  # Extract and return the CAT frequencies and rates
  
  # Open the iqtree file (if it exists)
  if (file.exists(iqtree_file) == TRUE){
    # Extract the best model
    best_model <- extract.best.model(iqtree_file)
    # Check if the best model is a CXX model
    cat_model_check <- grepl("C10|C20|C30|C40|C50|C60", best_model)
    
    # If the best model is not a CAT model, return NA
    # Otherwise, extract and return the frequency vector
    if (cat_model_check == FALSE){
      output_vector <- NA
    } else {
      # Extract and return the CXX frequency vectors
      # Open the .iqtree file:
      iq_lines <- readLines(iqtree_file)
      # Identify lines for table to extract
      table_start <- table(c(grep("No", iq_lines), grep("Component", iq_lines), grep("Rate", iq_lines), grep("Weight", iq_lines), grep("Parameters", iq_lines)))
      table_start_row <- as.numeric(names(which(table_start == 5)))
      empty_rows <- which(iq_lines == "")
      table_end_row <- empty_rows[table_start_row < empty_rows][1] - 1
      # Extract table of CXX parameters
      table_lines <- iq_lines[table_start_row:table_end_row]
      # Translate table into a table
      split_lines <- strsplit(table_lines, "  ")
      # Remove empty lines
      split_lines <- lapply(split_lines, function(x){x[x != ""]})
      # Turn strings into a dataframe
      cxx_table <- as.data.frame(do.call(rbind, split_lines[2:length(split_lines)]))
      names(cxx_table) <- gsub(" ", "", split_lines[[1]])
      cxx_table$Rate <- as.numeric(cxx_table$Rate)
      cxx_table$Weight <- as.numeric(cxx_table$Weight)
      # Check rates: output sum of weights
      print(paste0("IQ-Tree file: ", iqtree_file))
      print(paste0("Weight sum raw: ", sum(cxx_table$Weight)))
      sum_weights <- sum(as.numeric(cxx_table$Weight))
      # Ensure sum of weights = 1
      missing_weights = 1-sum_weights
      zero_weights = which(cxx_table$Weight == 0)
      if (length(zero_weights) > 0){
        # Replace missing weights (i.e. weight = 0) so sum of weights = 1
        fix_weights = missing_weights/length(zero_weights)
        cxx_table$Weight[zero_weights] <- fix_weights
      } else {
        # Add missing weights to biggest weight
        # Used to instead replace missing weights (i.e. weight = 0) so sum of weights = 1 (i.e if sum of weights = 0.95 and one component has weight 0, that component would be set to 0.05)
        # However, 0 weights are allowed in mixture models: see http://www.iqtree.org/release/v1.5.0/
        # Instead of adding rounding error to 0 weights, add it to the biggest weight (proportionally won't make as much of a difference)
        biggest_weight_row <- which(cxx_table$Weight == max(cxx_table$Weight))
        new_biggest_weight <- cxx_table$Weight[biggest_weight_row] + missing_weights
        cxx_table$Weight[biggest_weight_row] <- new_biggest_weight
      }
      # Reformat weights for nice output
      cxx_table$Weight <- format(as.numeric(cxx_table$Weight), digits = 4, scientific = FALSE)
      # Check rates: output sum of weights
      print(paste0("Weight sum post-check: ", sum(as.numeric(cxx_table$Weight))))
      # Paste the CXX parameters together into a single vector (to input into the MAST run)
      cxx_components <- paste0(cxx_table$Component, ":", cxx_table$Rate, ":", cxx_table$Weight)
      cxx_components <- gsub(" ", "", cxx_components)
      # Assemble into a CXX model
      cxx_model <- paste0("MIX{", paste(cxx_components, collapse = ","), "}")
      # Return the cxx model as the output_vector
      output_vector <- cxx_model
    } # end if (cat_model_check == FALSE)
  } else {
    output_vector <- NA
  } # end if (file.exists(iqtree_file) == TRUE)
  
  # Return output
  return(output_vector)
}



extract.alisim.model <- function(log_file){
  # Given a .log file (output from IQ-Tree), this function will extract the model specification from the Alisim instructions
  
  # Check if the log file exists
  if (file.exists(log_file) == TRUE){
    # If the log file does exist:
    ## Open the .log file:
    log_lines <- readLines(log_file)
    
    # Check for an Alisim line
    alisim_ind <- grep("ALISIM COMMAND", log_lines)
    # Check if the alisim_ind is present (if yes, that means there's a line containing Alisim simulation details)
    if (identical(alisim_ind,integer(0)) == FALSE & class(alisim_ind) == "integer"){
      ## Alisim model extraction
      alisim_line <- log_lines[grep("--alisim", log_lines)]
      # Split the line at the '\"' to get only the model 
      split_line <- strsplit(alisim_line, '\"')[[1]]
      # Need to only extract the part inside the '\"' - take the part of the split_line directly following "-m"
      #   e.g. if "-m" is in the first part of the line, take the second part of the line
      #  Because we split the line at the quotes, the part after "-m" will always be the model (as the model is in quotation marks)
      model_raw <- split_line[grep("-m", split_line)+1]
      # Remove any blank spaces from the mode
      alisim_model <- gsub(" ", "", model_raw)
    } else {
      # The Alisim line is missing, and we cannot report the model
      # Return NA for this file
      alisim_model <- NA
    } # end if (identical(rate_ind,logical(0)) == FALSE & class(rate_ind) == "integer")
  } # end if (file.exists(iqtree_file) == TRUE)
  
  # Return the output
  return(alisim_model)
}


extract.model.details <- function(iqtree_file){
  # Given a .iqtree file, this function will extract the model of sequence evolution parameters
  
  # read in the IQ-TREE file to get substitution model and parameters
  iq_file <- readLines(iqtree_file)
  # extract the file name
  ind      <- grep("Input file name:", iq_file)
  op1      <- gsub(" ", "", strsplit(iq_file[[ind]], ":")[[1]][2])
  # extract the number of taxa and extract the length of the alignment
  ind         <- grep("Input data:", iq_file)
  input_str   <- iq_file[[ind]] # get the line that contains this info
  input_ls    <- strsplit(input_str, " ")
  if (grepl("partitions", input_str) == FALSE){
    op2         <- input_ls[[1]][3] # extract number of sequences (number of taxa)
    op2.5       <- NA # Specify number of partitions
    op3         <- input_ls[[1]][6] # extract number of sites 
  } else if (grepl("partitions", input_str) == TRUE){
    op2         <- input_ls[[1]][3] # extract number of sequences (number of taxa)
    op2.5       <- input_ls[[1]][6] # Specify number of partitions
    op3         <- input_ls[[1]][9] # extract number of sites 
  }
  # Extract the model of substitution (same for amino-acid and nucleotide files)
  op4 <- extract.best.model(iqtree_file)
  # Extract information about the sequence alignment
  if (is.na(op2.5) == TRUE){
    # If the iqtree file is for a single alignment, extract information about the different kinds of sites
    ind <- grep("Number of constant sites:", iq_file)
    num_lines <- iq_file[c(ind:(ind + 3))] # take the four lines listing the number of different kinds of sites
    num_split <- unlist(strsplit(num_lines, ":")) # Split the lines at the colon
    num_names <- num_split[c(TRUE, FALSE)] # before the colon = name
    num_vals <- num_split[c(FALSE, TRUE)] # after the colon = value
    num_vals_regx <- regmatches(num_vals, gregexpr("[[:digit:]]+", num_vals)) # extract all the numbers after the colon
    # four strings = four lists of numbers: take first value from each list to get number of sites
    num_vals <- c(num_vals_regx[[1]][1], num_vals_regx[[2]][1], num_vals_regx[[3]][1], num_vals_regx[[4]][1]) 
  } else {
    # If the iqtree file is for multiple alignments, extract information about the different kinds of sites
    # Extract the lines containing information about the table of alignments
    table_ind <- grep("ID\tType\tSeq\tSite\tUnique\tInfor\tInvar\tConst\tName", iq_file) + 1
    end_table_ind <- grep("Column meanings:", iq_file) - 2
    # Extract rows in the table
    table_rows <- iq_file[table_ind:end_table_ind]
    # Split the table rows at the tabs ("\t")
    table_rows_split <- strsplit(table_rows, "\t")
    # Extract number of sites from each row ("Site")
    num_sites <- unlist(lapply(table_rows_split, "[[", 4))
    # Extract number of unique site patterns from each row ("Unique")
    num_unique <- unlist(lapply(table_rows_split, "[[", 5))
    # Extract number of parsimony-informative sites from each row ("Infor")
    num_infor <- unlist(lapply(table_rows_split, "[[", 6))
    # Extract number of invariant sites from each row ("Invar")
    num_invar <- unlist(lapply(table_rows_split, "[[", 7))
    # Extract number of constant sites from each row ("Const") - NOTE: can be subset of invariant sites
    num_const <- unlist(lapply(table_rows_split, "[[", 8))
    # Create the output vectors
    num_names <- c("Number of constant sites", "Number of invariant (constant or ambiguous constant) sites",
                   "Number of parsimony informative sites", "Number of distinct site patterns")
    num_vals <- c(sum(as.numeric(num_const)), sum(as.numeric(num_invar)), 
                  sum(as.numeric(num_infor)), sum(as.numeric(num_unique)) )
  } # end if (is.na(op2.5) == TRUE)
  
  # Check whether the sites are "amino-acid" or "nucleotide"
  # First, check the input_str to see if it specifies the alignment type
  if (grepl("amino-acid", input_str) == TRUE){
    input_type = "amino-acid"
  } else if (grepl("nucleotide", input_str) == TRUE){
    input_type = "nucleotide"
  } else {
    # If the alignment type is not specified in the input_str, check whether this is included in the partition table
    table_ind <- grep("ID\tType\tSeq\tSite\tUnique\tInfor\tInvar\tConst\tName", iq_file) + 1
    end_table_ind <- grep("Column meanings:", iq_file) - 2
    # Extract rows in the table
    table_rows <- iq_file[table_ind:end_table_ind]
    # Split the table rows at the tabs ("\t")
    table_rows_split <- strsplit(table_rows, "\t")
    # Extract the second element of each table row
    table_types <- unlist(lapply(table_rows_split, "[[", 2))
    # Check what the table_types are
    if (unique(table_types) == "AA"){
      # If all sequences are AA, then the input_type is AA
      input_type = "amino-acid"
    } else if (unique(table_types) == "DNA"){
      # If all sequences are DNA, then the input_type is DNA
      input_type = "nucleotide"
    } else {
      # If there is more than one type of sequence, then the sequences are a mix
      input_type = "mix"
    }
  }
  
  # if the input is DNA (nucleotide sites), gather that information
  if (input_type == "nucleotide"){
    # Check for presence of rate parameter details
    ind <- grep("A-C", iq_file)
    # Extract the rate parameters
    if (identical(ind, integer(0)) == FALSE){
      # Rate parameters are present in the file. Extract them
      rate1 <- as.numeric(strsplit(iq_file[[grep("A-C", iq_file)]],":")[[1]][2]) # A-C rate (same as code above, but combined 4 lines into 1 line)
      rate2 <- as.numeric(strsplit(iq_file[[grep("A-G", iq_file)]],":")[[1]][2]) # A-G rate
      rate3 <- as.numeric(strsplit(iq_file[[grep("A-T", iq_file)]],":")[[1]][2]) # A-T rate
      rate4 <- as.numeric(strsplit(iq_file[[grep("C-G", iq_file)]],":")[[1]][2]) # C-G rate
      rate5 <- as.numeric(strsplit(iq_file[[grep("C-T", iq_file)]],":")[[1]][2]) # C-T rate
      rate6 <- as.numeric(strsplit(iq_file[[grep("G-T", iq_file)]],":")[[1]][2]) # G-T rate
      # Create output vectors
      rate_names <- c("A-C_rate", "A-G_rate", "A-T_rate", "C-G_rate", "C-T_rate", "G-T_rate")
      rate_vals <- c(rate1, rate2, rate3, rate4, rate5, rate6)
    } else if (identical(ind, integer(0)) == TRUE){
      # Rate parameters are not present in the file
      # Create empty output vectors
      rate_names <- c()
      rate_vals <- c()
    }
    
    # Check for presence of state frequency details
    ind <- grep("State frequencies", iq_file)
    # Extract the state frequencies
    if (identical(ind, integer(0)) == FALSE){
      # State frequencies are present in the output file
      # Extract line of state frequencies
      state_freq_line <- iq_file[[ind]]
      if (state_freq_line == "State frequencies: (equal frequencies)"){
        # If the state frequencies are all equal, assign them all to 0.25 (1/4)
        sf1 <- 0.25 # pi(A) - A freq.
        sf2 <- 0.25 # pi(C) - C freq.
        sf3 <- 0.25 # pi(G) - G freq.
        sf4 <- 0.25 # pi(T) - T freq.
      } else {
        # If the state frequencies are not all equal, extract what they are
        sf1 <- as.numeric(strsplit(iq_file[[grep("pi\\(A\\)",iq_file)]], "=")[[1]][2]) # pi(A) - A freq. Remember to double backslash to escape before brackets
        sf2 <- as.numeric(strsplit(iq_file[[grep("pi\\(C\\)",iq_file)]], "=")[[1]][2]) # pi(C) - C freq.
        sf3 <- as.numeric(strsplit(iq_file[[grep("pi\\(G\\)",iq_file)]], "=")[[1]][2]) # pi(G) - G freq.
        sf4 <- as.numeric(strsplit(iq_file[[grep("pi\\(T\\)",iq_file)]], "=")[[1]][2]) # pi(T) - T freq.
      }
      # Create output vectors
      sf_names <- c("A_freq", "C_freq", "G_freq", "T_freq")
      sf_vals <- c(sf1, sf2, sf3, sf4)
    } else if (identical(ind, integer(0)) == TRUE){
      # State frequencies are present in the output file
      # Create output vectors
      sf_names <- c()
      sf_vals <- c()
    }
    
    # Check for presence of model of rate heterogeneity details
    ind <- grep("Model of rate heterogeneity:", iq_file)
    # Extract model of rate heterogeneity
    if (identical(ind, integer(0)) == FALSE){
      # If there is a section for rate heterogeneity, extract and output details 
      mrh1      <- strsplit(iq_file[[ind]], ":")[[1]][2] # Extract model of rate heterogeneity 
      mrh2      <- strsplit(iq_file[[ind+1]], ":")[[1]][2] # Line after the "model of rate heterogeneity" varies - extract it regardless of what it is 
      mrh2_name <- strsplit(iq_file[[ind+1]], ":")[[1]][1] # As the line varies, extract the name for the output dataframe
      mrh2_name <- gsub(" ", "_", mrh2_name) # change the name to be easy to parse
      # Assemble into output vector
      mrh_names <- c("model_of_rate_heterogeneity", "model_of_rate_heterogeneity_line2_name", "model_of_rate_heterogeneity_line2_value")
      mrh_vals <- c(mrh1, mrh2_name, mrh2)
    } else {
      # If there is no section for rate heterogeneity, ignore this section for output
      mrh_names <- c()
      mrh_vals <- c()
    }
    
    # make a list of the rows for the output dataframe
    names <- c("file_name", "sequence_type", "n_taxa", "npartitions", "n_sites", num_names, "substitution_model", 
               rate_names, sf_names, mrh_names)
    # Make a list of the output rows for the output dataframe
    op <- c(op1, "DNA", op2, op2.5, op3, num_vals, op4,
            rate_vals, sf_vals, mrh_vals)
    # Create the output dataframe
    op_df <- data.frame(names, op, stringsAsFactors = FALSE)
    # Name the columns
    names(op_df) <- c("parameter", "value")
    
    # Create the rate matrix Q
    Q_ind <- grep("Rate matrix Q:", iq_file)
    if (identical(Q_ind, integer(0)) == FALSE){
      # If the Q_ind exists, then the Q_df exists!
      # Extract the Q_df
      Q_start <- Q_ind + 2
      Q_end   <- Q_start + 3
      # Create the columns
      c1 <- c("A","C","G","T")
      c2 <- c()
      c3 <- c()
      c4 <- c()
      c5 <- c()
      # For each row in the iqtree file rate matrix
      for (i in Q_start:Q_end){
        # Split the row
        row <- strsplit(iq_file[[i]], " ")[[1]]
        row <- row[str_detect(row, "([0-9])")] # take only the numeric elements of the vector
        # Add the resulting values to the relevant columns
        c2 <- c(c2, as.numeric(row[1])) # convert to numeric so can use the numbers more easily later
        c3 <- c(c3, as.numeric(row[2]))
        c4 <- c(c4, as.numeric(row[3]))
        c5 <- c(c5, as.numeric(row[4]))
      }
      # Create a dataframe of the rate matrix Q
      q_df <- data.frame(c1, c2, c3, c4, c5, stringsAsFactors = FALSE)
      #Rename the columns
      names(q_df) <- c("nucleotide", "A", "C" ,"G", "T")
    } else {
      # If Q_ind is not present in this file, create an empty Q df
      q_df <- NA
    }
    
    # Check if the model of rate heterogeneity is present
    ind <- grep("Model of rate heterogeneity:", iq_file)
    if (identical(ind, integer(0)) == FALSE){
      # The model of rate heterogeneity is present
      # Check if the model for rate heterogeneity is uniform
      mrh1_check <- gsub(" ", "", mrh1)
      if (mrh1_check == "Uniform"){
        # If the model for rate heterogeneity is uniform, don't need to create a matrix for discrete gamma rate categories
        g_df <- "Uniform"
      } else {
        # If the model isn't uniform, need to create a matrix to collect and store the gamma category information
        #Create the matrix for discrete gamma categories
        g_start <- grep(" Category", iq_file) + 1 # get the index for the first line of the gamma categories matrix
        empty   <- which(iq_file == "") # get indexes of all empty lines
        empty   <- empty[empty > g_start] # get empty lines above gamma categories matrix
        g_end   <- empty[1] - 1 # get end index for gamma categories matrix (one less than next empty line)
        end_line <- iq_file[g_end]
        # if the end isn't an empty line, subtract one from the end count 
        # to exclude lines like "Relative rates are computed as MEAN of the portion of the Gamma distribution falling in the category."
        # to see if this is what's happening, check whether the line starts with a numeric section (i.e. a category for the gamma rate)
        check_line <- length(strsplit(strsplit(end_line, "        " )[[1]][1], " ")[[1]])
        if (check_line > 3){
          # If the check_line is longer than 3 characters, it won't be a group for the gamma categories but an instruction
          # Instructions can be excluded from the gamma matrix (but categories can't)
          g_end = g_end - 1
        }
        # Start collecting info for the matrix
        g1 <- c() # initialise columns to store data in
        g2 <- c()
        g3 <- c()
        # Iterate through rows in gamma matrix
        for (i in g_start:g_end){
          row <- strsplit(iq_file[[i]], " ")[[1]] # split the rows on the (large amount of) " "'s in the middle
          row <- row[str_detect(row, "([0-9])")] # take only the numeric elements of the vector
          g1 <- c(g1,as.numeric(row[1])) # add the values to the columns
          g2 <- c(g2,as.numeric(row[2]))
          g3 <- c(g3,as.numeric(row[3]))
        }
        g_df <- data.frame(g1, g2, g3, stringsAsFactors = FALSE) # create a dataframe of the information
        names(g_df) <- c("category", "relative_rate", "proportion") # name the columns
      }
    } else if (identical(ind, integer(0)) == TRUE){
      # The model of rate heterogeneity is not present. Return an empty dataframe.
      g_df <- NA
    }
    
    # Create a list of the three dataframes
    # This will be the output 
    params <- list(op_df, g_df, q_df)
    # Name the parameters so they're easy to access once you've outputted the data
    names(params) <- c("parameters", "gamma_categories", "Q_rate_matrix")
    
  } else if (input_type == "amino-acid"){ # alternatively if the data is amino acid sites
    # Check for presence of model of rate heterogeneity details
    ind <- grep("Model of rate heterogeneity:", iq_file)
    # Extract model of rate heterogeneity
    if (identical(ind, integer(0)) == FALSE){
      # If there is a section for rate heterogeneity, extract and output details 
      mrh1      <- strsplit(iq_file[[ind]], ":")[[1]][2] # Extract model of rate heterogeneity 
      mrh2      <- strsplit(iq_file[[ind+1]], ":")[[1]][2] # Line after the "model of rate heterogeneity" varies - extract it regardless of what it is 
      mrh2_name <- strsplit(iq_file[[ind+1]], ":")[[1]][1] # As the line varies, extract the name for the output dataframe
      # Assemble into output vector
      mrh_names <- c("model_of_rate_heterogeneity", "model_of_rate_heterogeneity_line2_name", "model_of_rate_heterogeneity_line2_value")
      mrh_vals <- c(mrh1, mrh2_name, mrh2)
    } else {
      # If there is no section for rate heterogeneity, ignore this section for output
      mrh_names <- c()
      mrh_vals <- c()
    }
    
    # Check for presence of state frequencies details
    ind <- grep("State frequencies:", iq_file)
    # Extract state frequencies
    if (identical(ind, integer(0)) == FALSE){
      # If there is a section for state frequencies, extract and output details 
      sf1      <- strsplit(iq_file[[ind]], ":")[[1]][2]
      # Assemble into output vector
      sf_names <- c("state_frequencies")
      sf_vals <- c(sf1)
    } else {
      # If there is no section for state frequencies, ignore this section for output
      sf_names <- c()
      sf_vals <- c()
    }
    
    # make a list of the rows for the output dataframe
    names <- c("file_name", "sequence_type", "n_taxa", "n_partitions", "n_sites", num_names, "substitution_model",
               mrh_names, sf_names)
    # Make a list of the output rows for the first output dataframe
    op <- c(op1,"amino-acid", op2, op2.5, op3, num_vals, op4, mrh_vals, sf_vals)
    op_df <- data.frame(names, op, stringsAsFactors = FALSE)
    names(op_df) <- c("parameter", "value")
    
    # Create a matrix of rate heterogeneity details (gamma categories):
    # Check for presence of model of rate heterogeneity details
    ind <- grep("Model of rate heterogeneity:", iq_file)
    if (identical(ind, integer(0)) == FALSE){
      # Check whether a gamma matrix is needed
      mrh1 <- strsplit(iq_file[[ind]], ":")[[1]][2] # Extract model of rate heterogeneity 
      mrh1_check <- gsub(" ", "", mrh1)
      if (mrh1_check == "Uniform"){
        # If the model for rate heterogeneity is uniform, don't need to create a matrix for discrete gamma rate categories
        g_df <- "Uniform"
      } else {
        #Create the matrix for discrete gamma categories
        g_start <- grep(" Category", iq_file) + 1 # get the index for the first line of the gamma categories matrix
        empty   <- which(iq_file == "") # get indexes of all empty lines
        empty   <- empty[empty > g_start] # get empty lines above gamma categories matrix
        g_end   <- empty[1] - 1 # get end index for gamma categories matrix (one less than next empty line)
        end_line <- iq_file[g_end]
        # if the end isn't an empty line, subtract one from the end count 
        # to exclude lines like "Relative rates are computed as MEAN of the portion of the Gamma distribution falling in the category."
        # to see if this is what's happening, check whether the line starts with a numeric section (i.e. a category for the gamma rate)
        check_line <- strsplit(strsplit(end_line, "        " )[[1]][1], " ")[[1]]
        check_line <- check_line[which(check_line != "")]
        check_line_count <- length(check_line)
        if (check_line_count > 3){
          # If the check_line is longer than 3 objects, it won't be a group for the gamma categories but an instruction
          # Eg. the relative rates line has 17 objects, because each word in the sentence is an object
          # Instructions can be excluded from the gamma matrix (but categories can't)
          g_end = g_end - 1
        }
        # initialise columns to store data in
        g1 <- c()
        g2 <- c()
        g3 <- c()
        # Iterate through rows in gamma matrix
        for (i in g_start:g_end){
          row <- strsplit(iq_file[[i]], "        ") # split the rows on the long strong of 0's in the middle
          g1 <- c(g1, as.numeric(row[[1]][1])) # add the values to the columns
          g2 <- c(g2, as.numeric(row[[1]][2]))
          g3 <- c(g3, as.numeric(row[[1]][3]))
        }
        g_df <- data.frame(g1, g2, g3, stringsAsFactors = FALSE) # create a dataframe of the information
        names(g_df) <- c("category", "relative_rate", "proportion") # name the columns 
      }
    } else if (identical(ind, integer(0)) == TRUE){
      # If no details on gamma categories are included, return nothing
      g_df <- NA
    }
    
    # Create a matrix of state frequency details:
    # Check for presence of state frequencies details
    ind <- grep("State frequencies:", iq_file)
    if (identical(ind, integer(0)) == FALSE){
      # If there is a section for state frequencies, extract and output details 
      sf1 <- strsplit(iq_file[[ind]], ":")[[1]][2]
      # Check whether state frequencies are needed
      sf1_squashed <- gsub(" ", "", sf1)
      if (sf1_squashed == "(empiricalcountsfromalignment)"){
        # Get starting line for frequencies
        start_ind <- grep("State frequencies:", iq_file) + 2
        # Take the 20 lines containing AA frequencies
        freq_lines <- iq_file[start_ind:(start_ind+19)]
        # Split up the frequency lines into the label and the frequency
        freq_split <- unlist(strsplit(freq_lines, "="))
        # Get the frequency
        freq_nums <- freq_split[c(FALSE, TRUE)]
        # Remove any spaces (from IQTree formatting)
        freq_nums <- gsub(" ","",freq_nums)
        # Get corresponding AA letter
        freq_names <- freq_split[c(TRUE, FALSE)]
        # Remove IQTree formatting
        freq_names <- gsub("pi\\(", "", freq_names)
        freq_names <- gsub("\\)", "", freq_names)
        freq_names <- gsub(" ", "", freq_names)
        # Create a nice dataframe
        f_df <- data.frame("amino_acid" = freq_names,
                           "frequency" = freq_nums,
                           stringsAsFactors = FALSE)
      } else {
        f_df <- "State frequencies from model"
      }
    } else if (identical(ind, integer(0)) == TRUE){
      # If no details on state frequencies are needed, return an empty dataframe
      f_df <- NA
    }
    
    # Create a list of the dataframes - this will be the output
    params <- list(op_df, g_df, f_df)
    # Name the parameters so they're easy to access once you've outputted the data
    names(params) <- c("parameters", "gamma_categories", "frequency")
  } 
  
  # Now the information has been collected, create an output dataframe
  # make the output dataframe
  return(params)
}


check.ModelFinder.models <- function(best_model, iqtree_file){
  # Function to take the best model by BIC for an alignment and check whether that model was tested by ModelFinder in IQ-Tree
  
  # Extract the mfp best model
  mfp_best_model <- extract.best.model(iqtree_file)
  # There's an output issue that sometimes means some models have "+I+I+", where it should be "+I+"
  #   Correct this in the best model and best model from ModelFinder, if applicable
  best_model <- gsub("\\+I\\+I\\+", "+I+", best_model)
  mfp_best_model <- gsub("\\+I\\+I\\+", "+I+", mfp_best_model)
  
  # Check if the iqtree file exists
  if (file.exists(iqtree_file) == TRUE){
    # If the iqtree_file does exist:
    ## Open the .iqtree file:
    iq_lines <- readLines(iqtree_file)
    
    ## Extract and format the ModelFinder section
    # Determine starting point of ModelFinder table
    start_line <- Reduce(intersect, list(grep("Model", iq_lines), grep("LogL", iq_lines), grep("AIC", iq_lines), grep("AICc", iq_lines), grep("BIC", iq_lines)))
    # Determine the end point of Modelfinder table
    # Find all the empty lines
    empty_lines <- which(iq_lines == "")
    end_line <- empty_lines[(which(empty_lines > start_line)[[1]])] - 1
    # Get the section containing ModelFinder lines
    mfp_lines <- iq_lines[start_line:end_line]
    # Split each line at the spaces
    mfp_lines_split <- strsplit(mfp_lines, " ")
    # Remove any empty characters, or any elements that are just "+" or "-"
    mfp_lines_split <- lapply(1:length(mfp_lines_split), function(i){mfp_lines_split[[i]][mfp_lines_split[[i]] != ""]})
    mfp_lines_split <- lapply(1:length(mfp_lines_split), function(i){mfp_lines_split[[i]][mfp_lines_split[[i]] != "-"]})
    mfp_lines_split <- lapply(1:length(mfp_lines_split), function(i){mfp_lines_split[[i]][mfp_lines_split[[i]] != "+"]})
    # Bind the second element of the list onwards into a nice dataframe
    mfp_df <- as.data.frame(do.call(rbind, mfp_lines_split[2:length(mfp_lines_split)]))
    # Use the first element of the list as the row names
    names(mfp_df) <- mfp_lines_split[[1]]
    # There's an output issue that sometimes means some models have "+I+I+", where it should be "+I+"
    #   Correct this issue in the mfp_df$Model column
    mfp_df$Model <- gsub("\\+I\\+I\\+", "+I+", mfp_df$Model)
    
    ## Investigate the models checked during model selection
    # Check whether the best_model was tested during model selection
    best_model_checked_bool <- (best_model %in% mfp_df$Model)
    # Get the model code from the best model (the first section of the model, e.g. "LG" for "LG+I+G4")
    best_model_code <- strsplit(best_model, "\\+")[[1]][1]
    # Check whether the model code from the best model was used at all 
    #     Reason for structure of logic: check whether any individual model in the dataframe contains that model_code 
    #     (i.e. if TRUE is present 1 or more times when we check "Is the model code in this model? True or false?")
    best_model_code_checked_bool <- (TRUE %in% (grepl(best_model_code, mfp_df$Model)))
    # Create the output vector
    mfp_check_op <- c(best_model_checked_bool, best_model_code_checked_bool)
    
  } else if (file.exists(iqtree_file) == FALSE){
    # If the iqtree_file doesn't exist, return NA for all entries in output vector
    mfp_check_op <- c(NA, NA)
  }
  
  # Return the output vector
  return(mfp_check_op)
}


extract.iqtree.file <- function(unique_id, IQTree_output_dir){
  # Small function to extract the relevant iqtree_file for a given alignment
  
  # Extract list of all files from the ML tree estimation 
  all_files <- list.files(IQTree_output_dir)
  # Extract the list of .iqtree files for this unique id
  iqtree_files <- grep(unique_id, grep("\\.iqtree", all_files, value = TRUE), value = TRUE)
  # Attach the full file path to the iqtree files
  iqtree_files <- paste0(IQTree_output_dir, iqtree_files)
  # Extract the MFP iqtree file
  mfp_iqtree_file <- grep("ModelFinder", iqtree_files, value = T)
  # Return the ModelFinder .iqtree file for this alignment
  return(mfp_iqtree_file)
}


extract.best.model.from.dataframe <- function(unique_id, best_models_df){
  # Small function to extract the relevant best_model for a given alignment
  
  # Extract dataset and matrix names from uid
  uid_dataset <- strsplit(unique_id, "\\.")[[1]][1]
  uid_matrix_name <- strsplit(unique_id, "\\.")[[1]][2]
  # Reduce dataframe to just rows containing the unique id
  uid_df <- best_models_df[best_models_df$dataset == uid_dataset & best_models_df$matrix_name == uid_matrix_name,]
  # Order matrix from highest to lowest BIC value
  uid_df <- uid_df[order(uid_df$best_model_BIC),]
  # Check whether there is an independent best model (2 rows) or whether MFP is the best model (1 row)
  if (nrow(uid_df) == 2){
    # Check whether the best model is from the ModelFinder run
    check_bool <- (uid_df$model_code[[1]] == "ModelFinder")
    # Create a vector of output values
    op_vec <- c(uid_dataset, uid_matrix_name, check_bool,
                uid_df$model_code[[1]], uid_df$best_model[[1]], uid_df$best_model_LogL[[1]], uid_df$best_model_BIC[[1]], uid_df$best_model_wBIC[[1]], 
                uid_df$tree_LogL[[1]], uid_df$tree_NumFreeParams[[1]], uid_df$tree_BIC[[1]], uid_df$tree_length[[1]], uid_df$tree_SumInternalBranch[[1]],
                uid_df$model_code[[2]], uid_df$best_model[[2]], uid_df$best_model_LogL[[2]], uid_df$best_model_BIC[[2]], uid_df$best_model_wBIC[[2]],
                uid_df$tree_LogL[[2]], uid_df$tree_NumFreeParams[[2]], uid_df$tree_BIC[[2]], uid_df$tree_length[[2]], uid_df$tree_SumInternalBranch[[2]])
  } else if (nrow(uid_df) == 1){
    # Check whether the best model is from the ModelFinder run
    check_bool <- (uid_df$model_code[[1]] == "ModelFinder")
    # Create a vector of output values
    op_vec <- c(uid_dataset, uid_matrix_name, check_bool,
                uid_df$model_code[[1]], uid_df$best_model[[1]], uid_df$best_model_LogL[[1]], uid_df$best_model_BIC[[1]], uid_df$best_model_wBIC[[1]], 
                uid_df$tree_LogL[[1]], uid_df$tree_NumFreeParams[[1]], uid_df$tree_BIC[[1]], uid_df$tree_length[[1]], uid_df$tree_SumInternalBranch[[1]],
                NA, NA, NA, NA, NA,
                NA, NA, NA, NA, NA)
  }
  # Return the output values
  return(op_vec)
}


check.ModelFinder.models.wrapper <- function(best_models_df, IQTree_output_dir){
  # Function to take in a single dataframe, and apply the check.ModelFinder.models function
  
  # Create the list of unique ids
  unique_ids <- unique(paste0(best_models_df$dataset, ".", best_models_df$matrix_name))
  
  # Extract information about best models, log likelihood, and BIC
  model_comparison_list <- lapply(unique_ids, extract.best.model.from.dataframe, best_models_df = best_models_df)
  # Convert list to dataframe
  model_comparison_df <- as.data.frame(do.call(rbind, model_comparison_list))
  names(model_comparison_df) <- c("dataset", "matrix_name", "best.model.is.ModelFinder",
                                  "best_model_code", "best_model_model", "best_model_LogL", "best_model_BIC", "best_model_wBIC",
                                  "best_model_tree_LogL", "best_model_tree_NumFreeParams", "best_model_tree_BIC", "best_model_tree_length", "best_model_tree_SumInternalBranch",
                                  "mfp_model_code", "mfp_model_model", "mfp_model_LogL", "mfp_model_BIC", "mfp_model_wBIC",
                                  "mfp_model_tree_LogL", "mfp_model_tree_NumFreeParams", "mfp_model_tree_BIC", "mfp_model_tree_length", "mfp_model_tree_SumInternalBranch")
  
  # Iterate through each dataset and extract the .iqtree file for the MFP run
  model_comparison_df$modelfinder_iqtree_files <- unlist(lapply(unique_ids, extract.iqtree.file, IQTree_output_dir = IQTree_output_dir))
  
  # Apply the function to one dataset at a time
  check_mfp_list <- lapply(1:nrow(model_comparison_df), function(i){check.ModelFinder.models(model_comparison_df$best_model_model[[i]], model_comparison_df$modelfinder_iqtree_files[[i]])})
  # Change output into a nince dataframe
  check_mfp_df <- as.data.frame(do.call(rbind,check_mfp_list))
  names(check_mfp_df) <- c("did.ModelFinder.test.best.model", "did.ModelFinder.test.best.model.code")
  # Bind the columns to the big dataframe
  model_comparison_df <- cbind(model_comparison_df, check_mfp_df)
  
  # Reorder the columns of the dataframe
  model_comparison_df <- model_comparison_df[c("dataset", "matrix_name", "best.model.is.ModelFinder","did.ModelFinder.test.best.model", "did.ModelFinder.test.best.model.code",
                                               "best_model_code", "best_model_model", "best_model_LogL", "best_model_BIC", "best_model_wBIC",
                                               "best_model_tree_LogL", "best_model_tree_NumFreeParams", "best_model_tree_BIC", "best_model_tree_length", "best_model_tree_SumInternalBranch",
                                               "mfp_model_code", "mfp_model_model", "mfp_model_LogL", "mfp_model_BIC", "mfp_model_wBIC",
                                               "mfp_model_tree_LogL", "mfp_model_tree_NumFreeParams", "mfp_model_tree_BIC", "mfp_model_tree_length", "mfp_model_tree_SumInternalBranch",
                                               "modelfinder_iqtree_files")]
  
  # Ensure all logical columns are logical class
  model_comparison_df$best.model.is.ModelFinder             <- as.logical(model_comparison_df$best.model.is.ModelFinder)
  model_comparison_df$did.ModelFinder.test.best.model       <- as.logical(model_comparison_df$did.ModelFinder.test.best.model)
  model_comparison_df$did.ModelFinder.test.best.model.code  <- as.logical(model_comparison_df$did.ModelFinder.test.best.model.code)
  
  # Return the dataframe of model information
  return(model_comparison_df)
}




#### Extract information from MAST runs ####
extract.HMM.output <- function(hmm_file){
  # Function to take an output prefix and directory, and return the results of the HMMster model
  
  # Create output label
  if (grepl("HMMster", hmm_file, ignore.case = T) == TRUE){
    MAST_option <- "HMMster"
  } else if (grepl("phyloHMM", hmm_file, ignore.case = T) == TRUE){
    MAST_option <- "phyloHMM"
  } else {
    MAST_option <- NA
  }
  # Open the iqtree file
  hmm_lines <- readLines(hmm_file)
  # Detect the output HMM probabilities
  hmm_prob_ind <- grep("Estimated HMM probabilities", hmm_lines, ignore.case = T)
  hmm_prob_line <- hmm_lines[ (hmm_prob_ind+1) ]
  hmm_probs <- strsplit(hmm_prob_line, "\t")[[1]]
  # Detect the number of sites for each category
  num_sites_ind <- grep("Number of sites for each category", hmm_lines, ignore.case = T)
  num_sites_line <- hmm_lines[ (num_sites_ind) ]
  num_sites_line_split <- strsplit(strsplit(num_sites_line, ":")[[1]][2], " ")[[1]]
  num_sites <- num_sites_line_split[which(num_sites_line_split != "")]
  # Detect the ratio of sites for each category
  rat_sites_ind <- grep("Ratio of sites for each category", hmm_lines, ignore.case = T)
  rat_sites_line <- hmm_lines[ (rat_sites_ind) ]
  rat_sites_line_split <- strsplit(strsplit(rat_sites_line, ":")[[1]][2], " ")[[1]]
  rat_sites <- rat_sites_line_split[which(rat_sites_line_split != "")]
  # Determine the number of trees
  num_trees <- length(hmm_probs)
  # Check how long each of the outputs are, and extend to 5 if necessary
  tws <- c(tws, rep(NA, (5 - num_trees )) )
  ttls <- c(ttls,rep(NA, (5 - num_trees )) )
  sibl <- c(sibl, rep(NA, (5 - num_trees )) )
  # Collect the output to return it
  hmm_output <- c(basename(hmm_file), MAST_option, num_trees, hmm_probs, num_sites, rat_sites)
  names(hmm_output) <- c("hmm_file", "analysis_type", "number_hypothesis_trees",paste0("tree_", 1:5, "_hmm_probs"),
                         paste0("tree_", 1:5, "_number_sites"), paste0("tree_", 1:5, "_ratio_sites"))
  # Return output
  return(hmm_output)
}


extract.tree.weights <- function(iq_file, trim.output.columns = FALSE){
  # Function to take an output prefix and directory, and return the results of the HMMster model
  
  # Open the iqtree file
  iq_lines <- readLines(iq_file)
  # Detect the tree weights for each tree
  tw_ind <- grep("Tree weights", iq_lines, ignore.case = T)
  tw_line <- iq_lines[ (tw_ind) ]
  tws <- gsub(" ", "", strsplit(strsplit(tw_line, ":")[[1]][2], ",")[[1]])
  # Detect the total tree length for each tree
  ttls_ind <- grep("Total tree lengths", iq_lines, ignore.case = T)
  ttls_line <- iq_lines[ (ttls_ind) ]
  ttls_raw <- gsub(" ", "", strsplit(strsplit(ttls_line, ":")[[1]][2], " ")[[1]])
  ttls <- ttls_raw[which(ttls_raw != "")]
  # Detect the sum of internal branch lengths for each tree
  sibl_ind <- grep("Sum of internal branch lengths", iq_lines, ignore.case = T)
  sibl_line <- iq_lines[ (sibl_ind) ]
  sibl_raw <- unlist(strsplit(strsplit(strsplit(sibl_line, ":")[[1]][2], "\\(")[[1]],  "\\)"))
  sibl <- gsub(" ", "", grep("\\%", sibl_raw, value = TRUE, invert = TRUE))
  # Determine the number of trees
  num_trees <- length(tws)
  # EITHER keep all 5 columns (for the maximum number of 5 trees) OR 
  #     remove any columns with NA values and return only the same number of columns as input trees
  if (trim.output.columns == FALSE){
    # Check how long each of the outputs are, and extend to 5 if necessary
    tws <- c(tws, rep(NA, (5 - num_trees )) )
    ttls <- c(ttls,rep(NA, (5 - num_trees )) )
    sibl <- c(sibl, rep(NA, (5 - num_trees )) )
    # Collect the output to return it
    mast_output <- c(basename(iq_file), num_trees, tws, ttls, sibl)
    names(mast_output) <- c("iq_file", "number_hypothesis_trees",paste0("tree_", 1:5, "_tree_weight"),
                            paste0("tree_", 1:5, "_total_tree_length"), paste0("tree_", 1:5, "_sum_internal_branch_lengths"))
  } else if (trim.output.columns == TRUE){
    mast_output <- c(basename(iq_file), num_trees, tws, ttls, sibl)
    names(mast_output) <- c("iq_file", "number_hypothesis_trees",paste0("tree_", 1:num_trees, "_tree_weight"),
                            paste0("tree_", 1:num_trees, "_total_tree_length"), paste0("tree_", 1:num_trees, "_sum_internal_branch_lengths"))
  }
  # Return output
  return(mast_output)
}




#### Extract information from maximum likelihood trees ####
# Functions to extract data from maximum likelihood trees estimated in IQ-Tree2 (from the .treefile output file)

check.tree.taxa <- function(tree_file){
  # Function to take a phylogenetic tree and return the list of taxa in that tree
  
  # Check whether the tree exists
  if (file.exists(tree_file) == TRUE){
    # Open the tree
    tree <- read.tree(file = tree_file)
    # Extract the list of taxa
    tree_tips <- tree$tip.label
  } else if (file.exists(tree_file) == FALSE){
    tree_tips <- NA
  }
  
  # Return the vector of tree tips
  return(tree_tips)
}


dataset.check.tree.taxa <- function(tree_files){
  # Function to take a list of tree files, check what taxa are in each tree, and determine whether the trees have the same taxa
  
  # Check how many tree files there are
  if (length(tree_files) > 1){
    # If there are multiple tree files, check if all the tree files have identical taxa
    
    # Trim any tree files that don't exist
    tree_files_exist <- tree_files[file.exists(tree_files)]
    
    # Apply the check.tree.taxa to each tree sequentially
    taxa_list <- lapply(tree_files_exist, check.tree.taxa)
    
    # Check whether all taxa have the same number of tips
    num_tree_tips <- unlist(lapply(taxa_list,length))
    num_taxa_lengths <- length(unique(num_tree_tips))
    
    # Check whether all trees have the same number of taxa
    if (num_taxa_lengths == 1){
      # All trees have the same number of taxa
      # Check whether the list of taxa is identical 
      # Sort each vector within the taxa_list so it is in alphabetical order
      taxa_list_alphabetical <- lapply(1:length(taxa_list), function(i){sort(taxa_list[[i]])})
      # Bind all tips into a dataframe (each row is a different tree)
      taxa_df <- as.data.frame(do.call(rbind, taxa_list_alphabetical))
      # Check whether all columns are identical
      #     FALSE = first originals; TRUE = duplicated 
      duplicated_check <- duplicated(taxa_df)
      # Count the number of TRUE and FALSE
      #   If all the `n` rows are identical, there should be 1 FALSE and (`n`-1) TRUE
      duplicated_check_table <- table(duplicated_check)
      num_original <- duplicated_check_table[["FALSE"]]
      num_copy <- duplicated_check_table[["TRUE"]]
      num_trees <- length(tree_files)
      
      # Check whether all trees have identical taxa
      if ((num_original == 1 & num_copy == (num_trees - 1)) | (num_original == 1 & num_copy == 0)){
        # There is one unique row (and all other rows are copies)
        # Therefore, all trees have identical taxa
        # Return the vector of taxa (by taking the first row of the dataframe as a character vector)
        ds_check_op <- as.character(taxa_df[1,])
      } else if (num_original > 1){
        # More than one one row is unique - different trees have different taxa
        # Return an error message
        ds_check_op <- "ERROR: different trees have different taxa"
      }
      
    } else if (num_taxa_lengths > 1){
      # Trees have different numbers of taxa      
      # Return an error message
      ds_check_op <- "ERROR: trees do not all have same number of taxa"
    }
  } else if (length(tree_files) == 1){
    # If there is a single tree file, return the taxa from that tree
    ds_check_op <- read.tree(file = tree_files)$tip.label
  }
  # Return the output (either a list of the taxa in all the trees, or an error message)
  return(ds_check_op)
}


dataset.check.tree.taxa.wrapper <- function(unique_ids, tree_folder){
  # Wrapper function for dataset.check.tree.taxa that operates on the completed maximum likelihood output dataframe
  
  # Each unique output id refers to an alignment: they are created by pasting the dataset and matrix name together (separated by a "." character)
  # Check that the unique ids are unique
  unique_ids <- unique(unique_ids)
  
  # Get the list of all ML trees from the maximim_likelihood_trees directory
  all_dir_files <- list.files(tree_folder)
  all_tree_files <- grep("treefile", all_dir_files, value = TRUE)
  # Paste the full directory path onto the basenames
  if (length(all_tree_files) > 0){
    all_tree_files <- paste0(tree_folder, all_tree_files)
  }
  
  # Iterate through the unique_ids and extract all tree files with that id
  dataset_tree_path_lists <- lapply(unique_ids, function(run_id){grep(run_id, all_tree_files, value = T)})
  
  # Determine if the trees for each unique_id have the same taxa
  dataset_taxa_lists <- lapply(1:length(dataset_tree_path_lists), function(i){dataset.check.tree.taxa(dataset_tree_path_lists[[i]])})
  
  # Change the output into a nice dataframe
  # Get the length of the longest element in the list 
  largest_length <- max(unlist(lapply(dataset_taxa_lists, length)))
  # Set lengths of all elements to the largest_length
  dataset_taxa_lists <- lapply(dataset_taxa_lists, function(v){ c(v, rep(NA, largest_length-length(v)))})
  # Bind list into a dataframe (each unique_id will be one column)
  dataset_taxa_df <- as.data.frame(do.call(cbind, dataset_taxa_lists))
  # Add names for each column
  names(dataset_taxa_df) <- unique_ids
  
  # Return the dataframe of names in each dataset/matrix combination
  return(dataset_taxa_df)
}




#### Extract information from tree mixture results ####
extract.tree.mixture.results <- function(tree_mixture_file, dataset, prefix, model, best_model, tree_branch_option){
  # Function to take one tree mixture file and output the results from the corresponding .iqtree file
  
  # Open tree mixture file
  iq_lines <- readLines(tree_mixture_file)
  
  # Extract input file details
  ind <- grep("Input file name", iq_lines)
  input_file_name <- gsub(" ", "", strsplit(iq_lines[ind], ":")[[1]][2])
  ind <- grep("User tree file name", iq_lines)
  input_tree_name <- gsub(" ", "", strsplit(iq_lines[ind], ":")[[1]][2])
  
  # Extract tree weights
  ind <- grep("Tree weights", iq_lines)
  tree_weight_str <- strsplit(iq_lines[ind], ":")[[1]][2]
  all_tree_weights <- strsplit(tree_weight_str, ",")[[1]]
  all_tree_weights <- as.numeric(gsub(" ", "", all_tree_weights))
  all_tree_weight_names <- paste0("tree_", 1:length(all_tree_weights), "_weight")
  
  # Extract tree lengths
  ind <- grep("Total tree lengths", iq_lines)
  tree_lengths_raw <- strsplit(strsplit(iq_lines[ind], ":")[[1]][2], " ")[[1]]
  tree_lengths <- tree_lengths_raw[tree_lengths_raw != ""]
  all_total_lengths <- paste0("tree_", 1:length(all_tree_weights), "_total_length")
  ind <- grep("Sum of internal branch lengths:", iq_lines)
  internal_lengths_raw <- strsplit(strsplit(iq_lines[ind], ":")[[1]][2], " ")[[1]]
  internal_lengths_trimmed <- grep("\\(|\\)|of|tree|length", internal_lengths_raw, invert = TRUE, value = TRUE)
  internal_lengths <- internal_lengths_trimmed[internal_lengths_trimmed != ""]
  all_internal_lengths <- paste0("tree_", 1:length(all_tree_weights), "_sum_internal_branch_lengths")
  internal_percents <- gsub("\\%", "", gsub("\\(", "", grep("\\%", internal_lengths_raw, value = TRUE)))
  all_internal_percent_names <- paste0("tree_", 1:length(all_tree_weights), "_internal_branch_percent")
  
  # Extract likelihood values
  ind <- grep("Log-likelihood of the tree", iq_lines)
  log_likelihood_str <- strsplit(iq_lines[ind], ":")[[1]][2]
  log_likelihood <- gsub(" ", "", strsplit(log_likelihood_str, "\\(")[[1]][1])
  log_likelihood_se <- gsub(" ", "", gsub(")", "", gsub("s.e.", "", strsplit(log_likelihood_str, "\\(")[[1]][2])))
  ind <- grep("Unconstrained log-likelihood", iq_lines)
  unconst_logl <- gsub(" ", "", strsplit(iq_lines[ind], ":")[[1]][2])
  ind <- grep("Number of free parameters", iq_lines)
  num_free_params <- gsub(" ", "", strsplit(iq_lines[ind], ":")[[1]][2])
  ind <- grep("Akaike information criterion", iq_lines)[1]
  aic_score <- gsub(" ", "", strsplit(iq_lines[ind], ":")[[1]][2])
  ind <- grep("Corrected Akaike information criterion", iq_lines)
  aicc_score <- gsub(" ", "", strsplit(iq_lines[ind], ":")[[1]][2])
  ind <- grep("Bayesian information criterion", iq_lines)
  bic_score <- gsub(" ", "", strsplit(iq_lines[ind], ":")[[1]][2])
  
  # Collate output vectors
  op_names <- c("prefix", "dataset", "model", "best_model", "tree_branch_option", 
                "log_likelihood", "log_likelihood_se", "unconstrained_log-likelihood",
                "number_free_parameters", "AIC_score", "AICc_score", "BIC_score",
                all_tree_weight_names, all_total_lengths, all_internal_lengths, all_internal_percent_names,
                "input_file_name", "input_tree_name", "tree_mixture_file")
  op_vals <- c(prefix, dataset, model, best_model, tree_branch_option, 
               log_likelihood, log_likelihood_se, unconst_logl, 
               num_free_params, aic_score, aicc_score, bic_score,
               all_tree_weights, tree_lengths, internal_lengths, internal_percents,
               input_file_name, input_tree_name, tree_mixture_file)
  # Assemble output dataframe
  op_df <- as.data.frame(matrix(op_vals, ncol = length(op_vals), nrow = 1, byrow = TRUE))
  names(op_df) <- op_names
  # Return output dataframe
  return(op_df)
}




#### Extract information from dataframes ####
check.remaining.runs <- function(dataset, input_parameter_file, output_parameter_file){
  ## Small function to check which (if any) runs are remaining for a single model
  # Open the input and output dataframes
  in_df <- read.table(input_parameter_file, header = TRUE, sep = "\t")
  out_df <- read.table(output_parameter_file, header = TRUE, sep = "\t")
  # Trim both dataframes to just the dataset of interest
  in_df <- in_df[grep(dataset, in_df$prefix),]
  out_df <- out_df[grep(dataset, out_df$prefix),]
  # Compare model codes
  in_model_codes <- in_df$model_code
  out_model_codes <- out_df$model_code
  # Check whether the model_code objects are identical
  if (setequal(in_model_codes, out_model_codes) == TRUE){
    # If the model_codes are identical, all runs are complete
    remaining_models = NA
    cat_only = NA
  } else {
    # If the model_codes are not identical, work out which models are missing
    missing_model_codes <- setdiff(in_model_codes, out_model_codes)
    # Collect missing model codes for output
    remaining_models <- paste(missing_model_codes, collapse = ",")
    # Check which of the remaining models are cat models
    detect_cat_models <- grep("C20|C60", grep("PMSF", missing_model_codes, value = T, invert = T), value = T)
    # Check whether all remaining models are CAT models
    if (identical(missing_model_codes, detect_cat_models)){
      # If all models are cat models, return cat_only = TRUE
      cat_only = TRUE
    } else {
      # If some/all remaining models are NOT cat models, return cat_only = FALSE
      cat_only = FALSE
    }
  }
  # Assemble and return output
  check_runs_op <- c(remaining_models, cat_only)
  names(check_runs_op) <- c("remaining_trees_to_run", "only_CXX_runs_remaining")
  return(check_runs_op)
}




#### Summarise information from output csvs ####
summarise.HMM.results <- function(row_index, hmm_df){
  ## Function to take one dataset and summarise the HMMster (MAST) results in one row
  # Get the relevant row from the dataset
  d_df <- hmm_df[row_index,]
  # Extract year of publication
  dataset_year <- as.numeric(str_extract(unique(d_df$dataset), "(\\d)+"))
  # Create a new row for the output
  new_row <- c(d_df$dataset, dataset_year, d_df$matrix_name, d_df$model_code,
               d_df$analysis_type, d_df$number_hypothesis_trees,
               d_df$tree_1_ratio_sites, d_df$tree_2_ratio_sites, 
               d_df$tree_3_ratio_sites, d_df$tree_4_ratio_sites, 
               d_df$tree_5_ratio_sites)
  names(new_row) <- c("dataset", "year", "matrix_name", "model_code",
                      "analysis_type", "number_hypothesis_trees", 
                      "tree_1_ratio_sites", "tree_2_ratio_sites",
                      "tree_3_ratio_sites", "tree_4_ratio_sites",
                      "tree_5_ratio_sites")
  # Return output
  return(new_row)
  
}

summarise.tree.weights <- function(row_index, tw_df){
  ## Function to take one dataset and summarise the HMMster (MAST) results in one row
  # Get the relevant row from the dataset
  d_df <- tw_df[row_index,]
  # Extract year of publication
  dataset_year <- as.numeric(str_extract(unique(d_df$dataset), "(\\d)+"))
  # Create a new row for the output
  new_row <- c(d_df$dataset, dataset_year, d_df$matrix_name, d_df$model_code,
               d_df$model_class, d_df$mast_branch_type,
               d_df$minimum_branch_length, d_df$number_hypothesis_trees,
               d_df$tree_1_tree_weight, d_df$tree_2_tree_weight, 
               d_df$tree_3_tree_weight, d_df$tree_4_tree_weight, 
               d_df$tree_5_tree_weight)
  names(new_row) <- c("dataset", "year", "matrix_name", "model_code",
                      "model_class", "mast_branch_type", 
                      "minimum_branch_length", "number_hypothesis_trees", 
                      "tree_1_tree_weight", "tree_2_tree_weight",
                      "tree_3_tree_weight", "tree_4_tree_weight",
                      "tree_5_tree_weight")
  # Return output
  return(new_row)
  
}

summarise.AU.test.results <- function(dataset_id, au_test_df){
  ## Function to take one dataset and summarise the AU test results in one row
  # Get rows of dataframe that have matching ID
  d_df <- au_test_df[au_test_df$ID == dataset_id,]
  # Extract year of publication
  dataset_year <- as.numeric(str_extract(unique(d_df$dataset), "(\\d)+"))
  # Identify the number of trees
  num_trees <- length(d_df$p_AU)
  # Pad the dataset df with NA (depending on the number of trees)
  au_test_pvalues <- c(d_df$p_AU, rep(NA, (5 - num_trees)) )
  # Create a new row for output
  new_row <- c(d_df$dataset[1], d_df$matrix[1], d_df$best_model_code[1], "AU_test_p_values", 
               au_test_pvalues, dataset_year)
  names(new_row) <- c("dataset", "matrix", "best_model_code", "topology_test", 
                      "tree_1", "tree_2", "tree_3", "tree_4", "tree_5", "year")
  # Return the output
  return(new_row)
}

summarise.eLW <- function(dataset_id, au_test_df){
  ## Function to take one dataset and summarise the AU test results in one row
  # Get rows of dataframe that have matching ID
  d_df <- au_test_df[au_test_df$ID == dataset_id,]
  # Extract year of publication
  dataset_year <- as.numeric(str_extract(unique(d_df$dataset), "(\\d)+"))
  # Identify the number of trees
  num_trees <- length(d_df$c_ELW)
  # Pad the dataset df with NA (depending on the number of trees)
  au_test_pvalues <- c(d_df$c_ELW, rep(NA, (5 - num_trees)) )
  # Create a new row for output
  new_row <- c(d_df$dataset[1], d_df$matrix[1], d_df$best_model_code[1], "expected_likelihood_weights", 
               au_test_pvalues, dataset_year)
  names(new_row) <- c("dataset", "matrix", "best_model_code", "topology_test", 
                      "tree_1", "tree_2", "tree_3", "tree_4", "tree_5", "year")
  # Return the output
  return(new_row)
}


summarise.topology.results <- function(dataset_id, topology_check_df, 
                                       excluded_models = c("C10", "C30", "C40", "C50")){
  ## Function to take one dataset and summarise the tree topology and Porifera clade topology
  # Get rows of dataframe that have matching ID
  d_df <- topology_check_df[grep(dataset_id, topology_check_df$ML_tree_name),]
  # Extract year of publication
  dataset_year <- as.numeric(str_extract(unique(d_df$dataset), "(\\d)+"))
  # Remove excluded models from the dataframe
  keep_rows <- which(!d_df$model_code %in% excluded_models)
  d_df <- d_df[keep_rows, ]
  # Remove any rows with NA or "" (i.e. missing trees)
  keep_tree_rows <- unique(c(which(d_df$sister_group != ""), which(d_df$PORI_topology != ""), 
                             which(!is.na(d_df$sister_group)), which(!is.na(d_df$PORI_topology))))
  d_df <- d_df[keep_tree_rows, ]
  # Summarise results
  topology_results <- c(length(which(d_df$sister_group == "Ctenophora")), 
                        length(which(d_df$sister_group == "Porifera")),
                        length(which(d_df$sister_group == "Ctenophora+Porifera")), 
                        length(which(d_df$sister_group == "Radiata")))
  pori_topology_results <- c(length(which(d_df$PORI_topology == "Monophyletic")), 
                             length(which(d_df$PORI_topology == "Paraphyletic")),
                             length(which(d_df$PORI_topology == "One taxon")))
  cten_cnid_mono_results <- c(length(which(d_df$`CTEN+CNID_monophyletic` == "Yes")), 
                              length(which(d_df$`CTEN+CNID_monophyletic` == "No")))
  # Summarise UFB support
  mean_UFB_support <- mean(as.numeric(d_df$UFB_support))
  # Check whether the dataset is completed
  num_trees <- nrow(d_df)
  num_completed_trees <- num_trees - length(which(d_df$sister_group == ""))
  if (num_completed_trees == 26){
    dataset_completed = TRUE
  } else {
    dataset_completed = FALSE
  }
  # Assemble output vector
  output_vector <- c(d_df$dataset[1], d_df$matrix_name[1], num_completed_trees, dataset_completed, topology_results, 
                     round((topology_results/num_trees*100), digits = 2), mean_UFB_support, pori_topology_results, 
                     round((pori_topology_results/num_trees*100), digits = 2), cten_cnid_mono_results, 
                     round((cten_cnid_mono_results/num_trees*100), digits = 2), d_df$PLAC_present[1], dataset_year)
  names(output_vector) <- c("dataset", "matrix_name", "num_completed_trees", "dataset_completed",
                            "CTEN_sister", "PORI_sister", "CTEN+PORI_sister", "Radiata_sister",
                            "percent_CTEN_sister", "percent_PORI_sister", "percent_CTEN+PORI_sister",
                            "percent_Radiata_sister", "mean_UFB_support", "PORI_monophyletic", 
                            "PORI_paraphyletic", "PORI_one_taxon", "percent_PORI_monophyletic",
                            "percent_PORI_paraphyletic", "percent_PORI_one_taxon", 
                            "CTEN+CNID_monophyletic", "CTEN+CNID_not_monophyletic",
                            "percent_CTEN+CNID_monophyletic", "percent_CTEN+CNID_not_monophyletic",
                            "PLAC_present", "dataset_year")
  # Return output vector
  return(output_vector)
}


tree.topology.results <- function(dataset_id, topology_check_df, model_order){
  ## Function to order the topologies consistently by model and return the results
  # Get rows of dataframe that have matching ID
  d_df <- topology_check_df[grep(dataset_id, topology_check_df$ML_tree_name),]
  # Remove excluded models from the dataframe
  keep_rows <- which(d_df$model_code %in% model_order)
  d_df <- d_df[keep_rows, ]
  # Sort by model order
  d_df <- d_df[match(model_order, d_df$model_code),]
  # Replace with shortened names
  d_df$sister_group_output <- d_df$sister_group
  d_df$sister_group_output[which(d_df$sister_group == "Ctenophora")] <- "CTEN"
  d_df$sister_group_output[which(d_df$sister_group == "Porifera")] <- "PORI"
  d_df$sister_group_output[which(d_df$sister_group == "Ctenophora+Porifera")] <- "CTEN+PORI"
  d_df$sister_group_output[which(d_df$sister_group == "Radiata")] <- "RADIATA"
  # Create output vector
  output_vector <- c(d_df$dataset[1], d_df$matrix_name[1], d_df$sister_group_output)
  # Return output
  return(output_vector)
}


porifera.topology.results <- function(dataset_id, topology_check_df, model_order){
  ## Function to order the topologies consistently by model and return the results
  # Get rows of dataframe that have matching ID
  d_df <- topology_check_df[grep(dataset_id, topology_check_df$ML_tree_name),]
  # Remove excluded models from the dataframe
  keep_rows <- which(d_df$model_code %in% model_order)
  d_df <- d_df[keep_rows, ]
  # Sort by model order
  d_df <- d_df[match(model_order, d_df$model_code),]
  # Replace with shortened names
  d_df$PORI_topology_output <- d_df$PORI_topology
  d_df$PORI_topology_output[which(d_df$PORI_topology == "Monophyletic")] <- "MONO"
  d_df$PORI_topology_output[which(d_df$PORI_topology == "Paraphyletic")] <- "PARA"
  d_df$PORI_topology_output[which(d_df$PORI_topology == "One taxon")] <- "NA"
  # Create output vector
  output_vector <- c(d_df$dataset[1], d_df$matrix_name[1], d_df$PORI_topology_output)
  # Return output
  return(output_vector)
}


return.MAST.branch.type <- function(file_path){
  # Quick function to return the branch type option from the MAST run
  
  if (grepl(".T.", file_path) == TRUE){
    branch_type = "T"
  } else if (grepl(".TR.", file_path) == TRUE){
    branch_type = "TR"
  } else {
    branch_type = NA
  }
  return(branch_type)
}



duplicate.constraint.rows <- function(df){
  ## Function to take in a dataframe and repeat each row a certain number of times  
  
  # Create the new dataframe using the first row
  first_row <- df[1, ]
  if (first_row$num_constraint_trees == 3){
    new_df <- rbind(first_row, first_row, first_row)
    new_df$constraint_tree_id <- 1:3
  } else if (first_row$num_constraint_trees == 5){
    new_df <- rbind(first_row, first_row, first_row, first_row, first_row)
    new_df$constraint_tree_id <- 1:5
  }
  # Iterate through the rows and duplicate the rows
  for (i in 2:nrow(df)){
    temp_row <- df[i, ]
    if (temp_row$num_constraint_trees == 3){
      temp_df <- rbind(temp_row, temp_row, temp_row)
      temp_df$constraint_tree_id <- 1:3
    } else if (temp_row$num_constraint_trees == 5){
      temp_df <- rbind(temp_row, temp_row, temp_row, temp_row, temp_row)
      temp_df$constraint_tree_id <- 1:5
    }
    new_df <- rbind(new_df, temp_df)
  }
  # Return the duplicated dataframe
  return(new_df)
}


