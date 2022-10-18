## caitlinch/metazoan-mixtures/func_data_processing.R
# Caitlin Cherryh 2022

# Functions for manipulating, processing, and preparing datasets

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


