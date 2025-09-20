# Set working directroy
setwd("C:/Users/aaisu/Box/Carter Lab/Projects/RNAseq_EthiopiaNUCI/Approach_II_outputs/Pathway_selection")

# Importing libraries
library(tidyverse)

# Creating a list to hold paths to each directory
directory_paths <- list()
directory_paths[["GO_enrichment"]] <- "C:/Users/aaisu/Box/Carter Lab/Projects/RNAseq_EthiopiaNUCI/Approach_II_outputs/GOenrichment/unfiltered_results"
directory_paths[["KEGGS"]] <- "C:/Users/aaisu/Box/Carter Lab/Projects/RNAseq_EthiopiaNUCI/Approach_II_outputs/KEGGS/unfiltered_results"
directory_paths[["GSEA"]] <- "C:/Users/aaisu/Box/Carter Lab/Projects/RNAseq_EthiopiaNUCI/Approach_II_outputs/GSEA/All_samples/Dataframes"

# Making a list to hold all unfiltered pathway and GO term data
all_data_nested <- list()

# Get the names from your directory_paths list. The loop will iterate through these.
directory_names <- names(directory_paths)

for (dir_name in directory_names) {
  # Get the actual file path corresponding to the current name
  current_path <- directory_paths[[dir_name]]
  
  cat(paste("Processing directory:", dir_name, "->", current_path, "\n"))
  
  # Find all files ending in .csv within the current directory.
  # We set `full.names = TRUE` to get the complete path for reading.
  csv_files <- list.files(path = current_path,
                          pattern = "\\.csv$",
                          full.names = TRUE)
  
  # Check if any CSV files were found in the directory.
  if (length(csv_files) == 0) {
    cat(paste("  -> No CSV files found in this directory. Skipping.\n"))
    # Use 'next' to skip the rest of the loop and move to the next directory
    next
  }
  
  # Read each CSV file into a temporary list of data frames.
  list_of_dataframes <- lapply(csv_files, function(file_path) {
    
    # Read the data from the file into a temporary data frame
    df <- read.csv(file_path)
    
    # Check if the 'p.adjust' column exists before trying to filter.
    # This makes the script robust and prevents errors if a file is missing the column.
    if ("p.adjust" %in% names(df)) {
      
      # Use subset() for a clean and readable filter.
      # This keeps only the rows where p.adjust is less than or equal to 0.05.
      # It also correctly handles any potential NA values in the column.
      df_filtered <- subset(df, p.adjust <= 0.05)
      
      # Return the newly filtered data frame
      return(df_filtered)
      
    } else {
      
      # If the column doesn't exist, print a warning and return the original, unfiltered data frame.
      cat(paste("  -> Warning: 'p.adjust' column not found in", basename(file_path), ". Returning original data.\n"))
      return(df)
      
    }
  })
  
  # Name each data frame in our temporary list using its original file name.
  # `basename()` is a handy function that strips the directory path, leaving just the file name.
  names(list_of_dataframes) <- basename(csv_files)
  
  # Assign the named list of data frames to our main list.
  # The key for this new entry is the name from your original `directory_paths` list.
  all_data_nested[[dir_name]] <- list_of_dataframes
  
  cat(paste("  -> Successfully imported", length(list_of_dataframes), "files.\n"))
}

# Creating a list to hold significantly up regulated and down regulated genes due to larval stages for each analysis
LS_pathways <- list()

# --- Extracting larval stage related pathways for each analysis and saving them in the list ---

# filtering out significantly up and down regulated enriched GO terms due to larval stages
go_enrichment_list <- all_data_nested$GO_enrichment
LS_pathways$GO_enrichment$LS_up <- as.data.frame(dplyr::bind_rows(all_data_nested$GO_enrichment[grepl("sigUp_IS2_vs_IS1|sigUp_IS3_vs_IS2|sigUp_IS4_vs_IS3", names(go_enrichment_list))]))
LS_pathways$GO_enrichment$LS_Down <- as.data.frame(dplyr::bind_rows(all_data_nested$GO_enrichment[grepl("sigDown_IS2_vs_IS1|sigDown_IS3_vs_IS2|sigDown_IS4_vs_IS3", names(go_enrichment_list))]))

# filtering out significantly up and down regulated enriched KEGGS pathways due to larval stages
KEGGS_list <- all_data_nested$KEGGS
LS_pathways$KEGGS$LS_up <- as.data.frame(dplyr::bind_rows(all_data_nested$KEGGS[grepl("Up_IS2_vs_IS1|Up_IS3_vs_IS2|Up_IS4_vs_IS3", names(KEGGS_list))]))
LS_pathways$KEGGS$LS_Down <- as.data.frame(dplyr::bind_rows(all_data_nested$KEGGS[grepl("Down_IS2_vs_IS1|Down_IS3_vs_IS2|Down_IS4_vs_IS3", names(KEGGS_list))]))

# filtering out significantly up and down regulated enriched GSEA pathways due to larval stages
GSEA_list <- all_data_nested$GSEA
LS_pathways$GSEA$LS_df <- as.data.frame(dplyr::bind_rows(all_data_nested$GSEA[grepl("IS2_vs_IS1|IS3_vs_IS2|IS4_vs_IS3", names(GSEA_list))]))

# --- Filtering out the larval stage related path ways from the other comparisons --- 

# Creating a nested data structure to hold the filtered pathways

filtered_pathways <- list()

# Running for loop to iterate over the unfiltered pathway data frames
for (analysis_name in names(all_data_nested)){
  
  # Implementing condition to select between GSEA and 
  if (analysis_name == "GSEA") {
    
    LS_df <- LS_pathways[[analysis_name]][["LS_df"]]
    
    for (df_name in names(all_data_nested[[analysis_name]])) {
      if (grepl("IS2_vs_IS1|IS3_vs_IS2|IS4_vs_IS3", df_name)) {next}
      filtered_IDs <- !(all_data_nested[[analysis_name]][[df_name]]$ID %in% LS_df$ID)
      filtered_pathways[[analysis_name]][[df_name]] <- all_data_nested[[analysis_name]][[df_name]][filtered_IDs,]
      if (nrow(filtered_pathways[[analysis_name]][[df_name]] > 0)) { 
        write.csv(filtered_pathways[[analysis_name]][[df_name]], file = paste("./", analysis_name, "/", df_name, sep = "")) 
      }
    }
    
  } else {
  
    LS_up_df <- LS_pathways[[analysis_name]][["LS_up"]]
    LS_down_df <- LS_pathways[[analysis_name]][["LS_Down"]]
    
    for (df_name in names(all_data_nested[[analysis_name]])) {
      # if (grepl("Up_IS2_vs_IS1|Up_IS3_vs_IS2|Up_IS4_vs_IS3|Down_IS2_vs_IS1|Down_IS3_vs_IS2|Down_IS4_vs_IS3", df_name)) {next}
      if (grepl("Up", df_name)){
        filtered_IDs <- !(all_data_nested[[analysis_name]][[df_name]]$ID %in% LS_up_df$ID)
        filtered_pathways[[analysis_name]][[df_name]] <- all_data_nested[[analysis_name]][[df_name]][filtered_IDs,]
      } else { 
        filtered_IDs <- !(all_data_nested[[analysis_name]][[df_name]]$ID %in% LS_down_df$ID)
        filtered_pathways[[analysis_name]][[df_name]] <- all_data_nested[[analysis_name]][[df_name]][filtered_IDs,]
      }
      if (nrow(filtered_pathways[[analysis_name]][[df_name]] > 0)) { 
        new_file_name <- sub("_unfiltered.csv", "", df_name)
        write.csv(filtered_pathways[[analysis_name]][[df_name]], file = paste("./", analysis_name, "/", new_file_name, ".csv", sep = "")) 
      }
    }
  }
}











