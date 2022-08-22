## caitlinch/metazoan-mixtures/func_constraint_trees.R
# Caitlin Cherryh 2022

# Functions for testing and applying constraint trees in iqtree2



#### Estimating an ML tree in IQ-Tree with a specified model ####
estimate.ml.iqtree <- function(iqtree2, alignment_file, model = "MFP", mset = NA, partition_file = NA, 
                               prefix = NA, number_parallel_threads = "AUTO", number_of_bootstraps = NA,
                               redo = FALSE, safe = FALSE, run.iqtree = TRUE){
  # Function to call iqtree and estimate a maximum likelihood tree using best practices
  
  # Add partition file if present
  if (is.na(partition_file) == TRUE){
    # If the partition file is NA, there is no partition file for this alignment
    partition_call <- ""
    # Check whether a model or mset is specified
    if (is.na(mset) == TRUE){
      # If mset = NA, then no mset option is specified.
      mset_call = ""
      # Tell IQ-Tree to use ModelFinder
      model_call = " -m MFP "
    } else if (is.na(mset) == FALSE){
      # If mset is specified, add mset command
      mset_call <- paste0(" -mset '", mset, "' ")
      # Do not use ModelFinder
      model_call = ""
    }
  } else if (is.na(partition_file) == FALSE){
    # If the partition file is not NA, add the command for a partition file to the command line for iqtree
    partition_call <- paste0(" -p ", partition_file, " ")
    # If there is a partition file, set the model selection to include a merging step
    model_call = " -m MFP+MERGE "
    # There is no mset command (models are already specified in the partition file)
    mset_call = ""
  } 
  
  # If prefix is specified, add a prefix command to the command line
  if (is.na(number_of_bootstraps) == FALSE){
    prefix_call = paste0(" -pre ", prefix, " ")
  } else if (is.na(number_of_bootstraps) == TRUE){
    prefix_call = ""
  }
  
  # If number of bootstraps is specified, add a bootstrap command to the command line
  if (is.na(number_of_bootstraps) == FALSE){
    bootstrap_call = paste0(" -B ", number_of_bootstraps, " ")
  } else if (is.na(number_of_bootstraps) == TRUE){
    bootstrap_call = ""
  }
  
  # If redo is TRUE, add redo command to IQ-Tree call
  if (redo == TRUE){
    redo_command = " -redo "
  } else if (redo == FALSE){
    redo_command = ""
  }
  
  # If safe is TRUE, add safe command to IQ-Tree call
  if (safe == TRUE){
    safe_command = " -safe "
  } else if (safe == FALSE){
    safe_command = ""
  }
  
  # Assemble the command line
  iqtree_call <- paste0(iqtree2, " -s ", alignment_file, partition_call, model_call, mset_call, prefix_call, 
                        " -nt ", number_parallel_threads, bootstrap_call, redo_command, safe_command)
  # Print the iqtree2 command
  print(iqtree_call)
  
  if (run.iqtree == TRUE){
    # Call iqtree to estimate the tree
    system(iqtree_call)
  } # end if (run.iqtree == TRUE)
  
} # end function



#### Extract details from IQ-Tree output files ####
extract.best.model <- function(iqtree_file){
  # Function that will extract the best model of sequence evolution or the model of sequence evolution used,
  #   given a .iqtree file
  
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
  
  # Return the best model (to be fed into further analyses)
  return(best_model)
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
  } else if (grepl("amino-acid", input_str) == TRUE){
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
      mrh2      <- as.numeric(strsplit(iq_file[[ind+1]], ":")[[1]][2]) # Line after the "model of rate heterogeneity" varies - extract it regardless of what it is 
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
        # If the model isn't uniform, need to create a matrix to collect and store the gamme category information
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
      mrh2      <- as.numeric(strsplit(iq_file[[ind+1]], ":")[[1]][2]) # Line after the "model of rate heterogeneity" varies - extract it regardless of what it is 
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
        check_line <- length(strsplit(strsplit(end_line, "        " )[[1]][1], " ")[[1]])
        if (check_line > 3){
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



#### Creating constraint trees ####
create.constraint.trees <- function(dataset, tree_id = NA, dataset_constraint_tree_dir, model, model_id, outgroup_taxa, ctenophora_taxa, 
                                    porifera_taxa, sponges_1_taxa, sponges_2_taxa, placozoa_taxa, cnidaria_taxa, bilateria_taxa,
                                    alignment_file, partitioned_check, partition_file, iqtree_path, number_parallel_threads){
  # Function to create the constraint trees and constraint tree information data frame, for a given dataset and model
  
  # Make sure you have an output id, which is a unique identifier for each dataset/alignment/model combination.
  if (is.na(tree_id) == FALSE){
    # If a tree_id is provided, use it in the file names
    output_id = tree_id
  }
  else if (is.na(tree_id) == TRUE){
    # If no tree_id is provided, create one
    output_id <- paste0(dataset, "_", model_id)
  }
  
  ## Hypothesis 1: Ctenophora-sister
  # Tree: (outgroup_taxa, (ctenophora_taxa, (porifera_taxa, (placozoa_taxa, cnidaria_taxa, bilateria_taxa))))
  # Construct constraint tree
  constraint_tree_1 <- paste0("((", 
                              paste(outgroup_taxa, collapse = ", "), 
                              "),((", 
                              paste(ctenophora_taxa, collapse = ", "), 
                              "),((", 
                              paste(c(porifera_taxa), collapse = ", "),
                              "), (",
                              paste(c(placozoa_taxa, cnidaria_taxa, bilateria_taxa), collapse = ", "),
                              "))));")
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, output_id, "_constraint_tree_", "1", ".nex")
  write(constraint_tree_1, file = constraint_tree_file_name)
  
  
  ## Hypothesis 2: Porifera-sister
  # Tree: (outgroup_taxa, (porifera_taxa, (ctenophora_taxa, (placozoa_taxa, cnidaria_taxa, bilateria_taxa))))
  # Construct constraint tree
  constraint_tree_2 <- paste0("((", 
                              paste(outgroup_taxa, collapse = ", "), 
                              "),((", 
                              paste(porifera_taxa, collapse = ", "), 
                              "),((", 
                              paste(c(ctenophora_taxa), collapse = ", "),
                              "),(",
                              paste(c(placozoa_taxa, cnidaria_taxa, bilateria_taxa), collapse = ", "), 
                              "))));")
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, output_id, "_constraint_tree_", "2", ".nex")
  write(constraint_tree_2, file = constraint_tree_file_name)
  
  
  ## Hypothesis 3: Porifera+Ctenophora-sister
  # Tree: (outgroup_taxa, ((porifera_taxa, ctenophora_taxa), (placozoa_taxa, (cnidaria_taxa, bilateria_taxa))))
  # Construct constraint tree
  constraint_tree_3 <- paste0("((", 
                              paste(outgroup_taxa, collapse = ", "), 
                              "),((", 
                              paste(c(porifera_taxa, ctenophora_taxa), collapse = ", "), 
                              "),(", 
                              paste(c(placozoa_taxa), collapse = ", "),
                              ", (",
                              paste(c(cnidaria_taxa, bilateria_taxa), collapse = ", "), 
                              "))));")
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, output_id, "_constraint_tree_", "3", ".nex")
  write(constraint_tree_3, file = constraint_tree_file_name)
  
  ## Hypothesis 4: Paraphyletic sponges, Porifera-sister
  # Tree: (outgroup_taxa, (sponges_1_taxa, (sponges_2_taxa, (ctenophora_taxa, (placozoa_taxa, cnidaria_taxa, bilateria_taxa)))))
  # Construct constraint tree
  constraint_tree_4 <- paste0("((", 
                              paste(outgroup_taxa, collapse = ", ") ,
                              ") ,((", 
                              paste(sponges_1_taxa, collapse = ", "), 
                              "), ((", 
                              paste(sponges_2_taxa, collapse = ", "), 
                              "), ((", 
                              paste(c(ctenophora_taxa), collapse = ", "),
                              "), (",
                              paste(c(placozoa_taxa, cnidaria_taxa, bilateria_taxa), collapse = ", "),
                              ")))));")
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, output_id, "_constraint_tree_", "4", ".nex")
  write(constraint_tree_4, file = constraint_tree_file_name)
  
  ## Hypothesis 5: Paraphyletic sponges, Ctenophora-sister
  # Tree: (outgroup_taxa, (ctenophora_taxa, (sponges_1_taxa, (sponges_2_taxa, (placozoa_taxa, cnidaria_taxa, bilateria_taxa)))))
  # Construct constraint tree
  constraint_tree_5 <- paste0("((", 
                              paste(outgroup_taxa, collapse = ", "),
                              ") ,((",
                              paste(ctenophora_taxa, collapse = ", "),
                              "), ((", 
                              paste(sponges_1_taxa, collapse = ", "),
                              "), ((", 
                              paste(c(sponges_2_taxa), collapse = ", "),
                              "), (",
                              paste(c(placozoa_taxa, cnidaria_taxa, bilateria_taxa), collapse = ", "),
                              ")))));")
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, output_id, "_constraint_tree_", "5", ".nex")
  write(constraint_tree_5, file = constraint_tree_file_name)
  
  ## Hypothesis 6: Paraphyletic sponges, Porifera-sister
  # Tree: (outgroup_taxa, (sponges_2_taxa, (sponges_1_taxa, (ctenophora_taxa, (placozoa_taxa, cnidaria_taxa, bilateria_taxa)))))
  # Construct constraint tree
  constraint_tree_6 <- paste0("((", 
                              paste(outgroup_taxa, collapse = ", ") ,
                              ") ,((", 
                              paste(sponges_2_taxa, collapse = ", "), 
                              "), ((", 
                              paste(sponges_1_taxa, collapse = ", "), 
                              "), ((", 
                              paste(c(ctenophora_taxa), collapse = ", "),
                              "), (",
                              paste(c(placozoa_taxa, cnidaria_taxa, bilateria_taxa), collapse = ", "),
                              ")))));")
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, output_id, "_constraint_tree_", "6", ".nex")
  write(constraint_tree_6, file = constraint_tree_file_name)
  
  ## Hypothesis 7: Paraphyletic sponges, Ctenophora-sister
  # Tree: (outgroup_taxa, (ctenophora_taxa, (sponges_2_taxa, (sponges_1_taxa, (placozoa_taxa, cnidaria_taxa, bilateria_taxa)))))
  # Construct constraint tree
  constraint_tree_7 <- paste0("((", 
                              paste(outgroup_taxa, collapse = ", "),
                              ") ,((",
                              paste(ctenophora_taxa, collapse = ", "),
                              "), ((", 
                              paste(sponges_2_taxa, collapse = ", "),
                              "), ((", 
                              paste(c(sponges_1_taxa), collapse = ", "),
                              "), (",
                              paste(c(placozoa_taxa, cnidaria_taxa, bilateria_taxa), collapse = ", "),
                              ")))));")
  constraint_tree_file_name <- paste0(dataset_constraint_tree_dir, output_id, "_constraint_tree_", "7", ".nex")
  write(constraint_tree_7, file = constraint_tree_file_name)
  
  # Assemble dataframe of information about the constraint trees
  constraint_df <- data.frame(constraint_tree_id = 1:7,
                              constraint_tree_paths = paste0(dataset_constraint_tree_dir, output_id, "_constraint_tree_", 1:7, ".nex"),
                              constraint_prefixes = paste0(output_id, "_ML_H", 1:7),
                              alignment_path = alignment_file,
                              model = model,
                              iqtree_path = iqtree_path,
                              constraint_trees = c(constraint_tree_1, constraint_tree_2, constraint_tree_3, 
                                                   constraint_tree_4, constraint_tree_5, constraint_tree_6,
                                                   constraint_tree_7),
                              num_threads = number_parallel_threads,
                              partitioned = partitioned_check,
                              partition_file = partition_file)
  
  # Write dataframe of information about constraint trees
  constraint_df_path <- paste0(dataset_constraint_tree_dir, output_id, "_constraint_tree_parameters.csv")
  write.csv(constraint_df, constraint_df_path, row.names = FALSE)
  
  # Return the constraint tree dataframe
  return(constraint_df)
}



#### Estimating ML trees using a constraint tree ####
run.iqtree.with.constraint.tree <- function(alignment_path, constraint_tree_file, partitioned_check = FALSE, partition_file = NA, 
                                            iqtree_path = "iqtree2", prefix = NA, model = NA, num_threads = 1, run.iqtree = TRUE){
  # Function to apply IQ-Tree to a series of alignments with a constraint tree
  
  # Set model for IQ-Tree run
  if (is.na(model) == TRUE){
    # If no model specified for IQ-Tree, use model finder (-m MFP) command
    iq_model = "MFP"
  } else {
    # Otherwise, use model specified in function call
    iq_model = model
  }
  
  # Add partition file if present
  if (partitioned_check == FALSE){
    partition_call <- ""
  } else if (partitioned_check == TRUE){
    partition_call <- paste0(" -p ", partition_file, " ")
  } 
  
  # Add prefix if present
  if (is.na(prefix) == TRUE){
    # If prefix is NA, do nothing
    prefix_call <- ""
  } else if (is.na(prefix) == FALSE){
    # If prefix is NA, add prefix to command line 
    prefix_call <- paste0(" --prefix ", prefix, " ")
  } 
  
  # Collate iqtree command
  iqtree_call <- paste0(iqtree_path, " -s ", alignment_path,  partition_call, " -m ", iq_model, " -g ", constraint_tree_file, " -T ", num_threads,  prefix_call)
  
  # Print IQ-Tree command
  print(iqtree_call)
  
  if (run.iqtree == TRUE){
    # Run IQ-Tree
    system(iqtree_call)
  } # end if (run.iqtree == TRUE)
  
} # end function

run.one.constraint.tree <- function(index, df, run.iqtree = TRUE){
  # Quick function to take in a dataframe, take relevant variables, and call the run.iqtree.with.constraint.tree function
  
  # Identify row
  row <- df[index, ]
  
  # Feed row information into function call
  # Call function with lapply whether run.iqtree = TRUE or run.iqtree = FALSE:
  #   either way, want to run function to print iqtree command line
  run.iqtree.with.constraint.tree(alignment_path = row$alignment_path, constraint_tree_file = row$constraint_tree_paths, 
                                  partitioned_check = row$partitioned_check, partition_file = row$partition_file, 
                                  iqtree_path = row$iqtree_path, prefix = row$constraint_prefixes, model = row$model,
                                  num_threads = row$num_threads, run.iqtree = run.iqtree)
}

run.one.constraint.dataframe <- function(csv_file, run.iqtree = TRUE){
  # Quick function to take in a dataframe, and estimate hypothesis trees by feeding it row by row into the run.one.constraint.tree function
  
  # Open the dataframe
  df <- read.csv(csv_file)
  # Estimate an ML tree in IQ-Tree for each constraint tree
  # Call function with lapply whether run.iqtree = TRUE or run.iqtree = FALSE:
  #   either way, want to run function to print iqtree command line
  lapply(1:nrow(df), run.one.constraint.tree, df, run.iqtree = run.iqtree)
}



#### Collating multiple trees into a single file ####
combine.hypothesis.trees <- function(tree_id, constraint_tree_directory, outgroup_taxa = NA){
  # Function to open all hypothesis trees with a given id in a folder and collate them into one file
  
  # List all hypothesis trees in the folder
  all_constraint_tree_dir_files <- list.files(constraint_tree_directory, recursive = TRUE)
  # Remove any files with "ignore" in the name
  all_constraint_tree_dir_files <- grep("ignore", all_constraint_tree_dir_files, value = TRUE, invert = TRUE)
  # Find all files for this tree_id
  tree_id_files <- grep(tree_id, all_constraint_tree_dir_files, value = TRUE)
  # Find all hypothesis trees for this tree_id (hypothesis trees are marked by HX, where 1<= X <= 7)
  hypothesis_tree_files <- grep("H1|H2|H3|H4|H5|H6|H7", tree_id_files, value = TRUE)
  # Extend file path
  if (length(hypothesis_tree_files) > 0){
    hypothesis_tree_files <- paste0(constraint_tree_directory, hypothesis_tree_files)
  }
  
  # Read in hypothesis tree files
  hypothesis_trees <- lapply(hypothesis_tree_files, read.tree)
  # Convert hypothesis_trees from a list into a multiPhylo object 
  class(hypothesis_trees) <- "multiPhylo"
  
  # Output the (unrooted) hypothesis trees
  unrooted_file <- paste0(constraint_tree_directory, tree_id, "_unrooted_hypothesis_trees.tre")
  write.tree(hypothesis_trees, file = unrooted_file)
  
  if (is.na(outgroup_taxa) == FALSE){
    # If the outgroup taxa are provided, root the hypothesis trees and save the rooted trees
    # Root hypothesis trees
    rooted_hypothesis_trees <- root(hypothesis_trees, outgroup_taxa)
    # Output the rooted hypothesis trees
    rooted_file <- paste0(constraint_tree_directory, tree_id, "_rooted_hypothesis_trees.tre")
    write.tree(hypothesis_trees, file = rooted_file)
    # Return paths for both rooted and unrooted hypothesis trees
    op_vec <- c(rooted_file, unrooted_file)
    names(op_vec) <- c("rooted_hypothesis_tree_file", "unrooted_hypothesis_tree_file")
  } else if (is.na(outgroup_taxa) == TRUE){
    # If the outgroup taxa are not provided, return only the path to the unrooted hypothesis trees
    op_vec <- c(unrooted_file)
    names(op_vec) <- c("unrooted_hypothesis_tree_file")
  }
  
  # Return the file names
  return(op_vec)
}



#### Applying the mixture of trees model ####
run.tree.mixture.model <- function(alignment_file, hypothesis_tree_file, partition_file, use.partition = FALSE, 
                                   prefix, model, iqtree2_tree_mixtures_implementation, tree_branch_option = "TR",
                                   number_parallel_threads, run.iqtree = TRUE){
  # Function runs the IQ-Tree2 mixture of trees model implementation given a sequence alignment, a set of hypothesis trees, and details about the model.
  # Currently cannot run with partition model
  # Tree branch options: 
  #   tree_branch_option = "TR" <- trees will have same length branches (for branches that appear in 2 or more trees)
  #   tree_branch_option = "T"  <- tree branches can have different lengths on different trees
  
  # Add model if present
  if (is.na(model) == TRUE){
    # If no model specified for IQ-Tree, use model finder (-m MFP) command
    if (is.na(partition_file) == TRUE){
      # If no partition file is present, use only the MFP command
      model_call = "MFP"
    } else if (is.na(partition_file) == FALSE){
      # If a partition file is present, use MFP+MERGE
      model_call = "MFP+MERGE"
    }
    # Extend the model to have the +TR command 
    model_call = paste0("'",model_call, "+TR'")
  } else if (is.na(model) == FALSE){
    # If a model is provided, use that model
    model_call = model
    # Extend the model to have the +TR command 
    model_call = paste0("'", model_call, "+", tree_branch_option, "'")
  }
  
  # Add partition file if present
  if (is.na(partition_file) == TRUE){
    # If partition_file is NA, do nothing
    partition_call <- ""
  } else if (is.na(partition_file) == FALSE){
    # If prefix is NA, add prefix to command line 
    partition_call <- paste0("-p ", partition_file,)
    # Pad partition_call with white space (for pasting into command line)
    partition_call <- paste0(" ", partition_call, " ")
  } 
  
  # Add prefix if present
  if (is.na(prefix) == TRUE){
    # If prefix is NA, do nothing
    prefix_call <- ""
  } else if (is.na(prefix) == FALSE){
    # If prefix is NA, add prefix to command line 
    prefix_call <- paste0("-pre ", prefix)
    # Pad prefix_call with white space (for pasting into command line)
    prefix_call <- paste0(" ", prefix_call, " ")
  }
  
  if (use.partition == FALSE){
    # Assemble the command for the tree mixtures model
    treemix_command <- paste0(iqtree2_tree_mixtures_implementation, " -s ", alignment_file, 
                              " -m  ", model_call, " -te ", hypothesis_tree_file, 
                              " -nt ", number_parallel_threads, prefix_call)
  } else if (use.partition == TRUE){
    # Assemble the command for the tree mixtures model
    treemix_command <- paste0(iqtree2_tree_mixtures_implementation, " -s ", alignment_file, partition_call, 
                              " -m ", model_call, " -te ", hypothesis_tree_file, 
                              " -nt ", number_parallel_threads, prefix_call)
  }
  
  # Change working directories (to store IQ-Tree output files in the right place)
  setwd(dirname(hypothesis_tree_file))
  
  # Print the iqtree command line
  print(treemix_command)
  
  if (run.iqtree == TRUE){
    # Call IQ-Tree2 with the command
    system(treemix_command)
  } # end if (run.iqtree == TRUE)
  
} # end function



#### Extract information from tree mixture results ####
extract.tree.mixture.results <- function(tree_mixture_file){
  # Function to take one tree mixture file and output the results from the corresponding .iqtree file
  
}



