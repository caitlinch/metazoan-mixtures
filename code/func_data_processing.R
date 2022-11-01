## caitlinch/metazoan-mixtures/func_data_processing.R
# Caitlin Cherryh 2022

# Functions for manipulating, processing, and preparing datasets

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
  }
  
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

