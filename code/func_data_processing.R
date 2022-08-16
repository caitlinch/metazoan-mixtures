## caitlinch/metazoan-mixtures/func_data_processing.R
# Caitlin Cherryh 2022

# Functions for manipulating, processing, and preparing datasets

extract.partition.models <- function(partition_file, paper_id){
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
  # Identify how many times each model was used
  model_t <- table(models)
  model_count <- as.numeric(model_t)
  model_name <- names(model_t)
  # Construct a dataframe with the model name and information
  model_df <- data.frame(paper_id = paper_id,
                         partition_file = basename(partition_file),
                         model = model_name,
                        count = model_count)
  # Get the starting names from the model_df
  model_df_init_names <- names(model_df)
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
  # Check whether each model contains +I, +F, or +G. Add results to model dataframe
  model_df$plus_f <- grepl("\\+F", model_name)
  model_df$plus_i <- grepl("\\+I\\+", model_name)
  model_df$plus_g <- grepl("\\+G", model_name)
  # Return model dataframe
  return(model_df)
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


