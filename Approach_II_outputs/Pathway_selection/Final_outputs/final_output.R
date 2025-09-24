# Set working directroy
setwd("C:/Users/aaisu/Box/Carter Lab/Projects/RNAseq_EthiopiaNUCI/Approach_II_outputs/Pathway_selection/Final_outputs")

# Importing libraries
library(tidyverse)
library(VennDiagram)

# Creating a list to hold paths to each directory containing CSV files of KEGGS, GO enrichment and GSEA results with unfiltered significant genes
directory_paths <- list()
directory_paths[["GO_enrichment"]] <- "../GO_enrichment/"
directory_paths[["KEGGS"]] <- "../KEGGS/"
directory_paths[["GSEA"]] <- "../GSEA/"

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

# --- Organizing the data for creating Venn diagrams ---

# --- 1. Define Your Combinations ---
# This list tells the script which data frames to combine.
# The combinations were saved in the CSV file with columns "union_name,analysis_name,df_names,output_name,id_column"
# You can add as many combinations as you need.
# The for loop is populating the list from the CSV file

# Reading in combinations
combinations_csv <- read.csv("combinations.csv", header = TRUE, sep = ",")

# Creating list to hold the combinations
combinations_to_create <- list()

# Populating the list
for (i in 1:nrow(combinations_csv)){
  # print(combinations_csv[i,3])
  combinations_to_create[[combinations_csv[i,1]]] <- list(
    analysis = combinations_csv[i,2],
    dfs = unlist(strsplit(combinations_csv[i,3], , split = "\\s*\\|\\s*")),
    output_name = combinations_csv[i,4],
    id_column = combinations_csv[i,5]
  )
}


# --- 2. Process Combinations and Create Unions ---
# This loop iterates through your definitions and performs the main task.

for (combo in combinations_to_create) {
  
  # Get the list of data frames for the current analysis type (e.g., all_data_nested$KEGGS)
  target_dfs_list <- all_data_nested[[combo$analysis]]
  
  # Select only the data frames specified in the current combination
  # `na.omit()` ensures we skip any names that weren't found in the list
  dfs_to_merge <- target_dfs_list[na.omit(match(combo$dfs, names(target_dfs_list)))]
  
  # Check if we have any data frames to merge before proceeding
  if (length(dfs_to_merge) > 0) {
    
    # Use bind_rows() to stack all selected data frames into one large data frame.
    # It intelligently handles cases where columns might not match.
    combined_df <- dplyr::bind_rows(dfs_to_merge)
    
    # Find the "union" by keeping only the unique rows based on your specified ID column.
    # `distinct()` keeps the first row it encounters for each unique ID.
    union_df <- combined_df %>%
      dplyr::distinct(!!sym(combo$id_column), .keep_all = TRUE)
    
    # Save the new union data frame back into the nested list
    all_data_nested[[combo$analysis]][[combo$output_name]] <- union_df
    
    cat(paste("Created union '", combo$output_name, "' with ", nrow(union_df), " unique entries.\n", sep = ""))
    
  } else {
    warning(paste("Could not find any of the specified data frames for '", combo$output_name, "'. Skipping.", sep = ""))
  }
}

# Clean up variables from the loop
rm(combo, combinations_to_create)

### --- Creating Venn Diagrams --- ###

# --- 1. Define the Venn Diagrams to Create ---
# This list is your configuration. Each entry defines one Venn diagram.
# For each diagram, specify:
#   - dfs_to_compare: A list where each entry points to a data frame.
#       - analysis: The top-level name (e.g., "KEGGS").
#       - df_name: The name of the data frame (e.g., "Union_LarvalStage_Up").
#   - title: The title that will appear on the Venn diagram chart.
#   - filename: The name of the image file to save (e.g., "my_venn.png").
#
# NOTE: The VennDiagram package works best for 2, 3, 4, or 5 sets.

venn_definitions <- list(
  
  # --- A 3-set Venn for Up-regulated Wild Caught Larvae KEGGS ---
  sites_keggs = list(
    dfs_to_compare = list(
      list(analysis = "KEGGS", df_name = "sigUp_Adama_U"),
      list(analysis = "KEGGS", df_name = "sigUp_Erer_U"),
      list(analysis = "KEGGS", df_name = "sigUp_Jijiga_U")
    ),
    title = "KEGG Pathways shared among sites",
    filename = "Venn_Sites_KEGGS.png"
  ),
  
  # --- A 3-set Venn for Up-regulated Wild Caught Larvae GO enrichments ---
  sites_GO_enrichment = list(
    dfs_to_compare = list(
      list(analysis = "GO_enrichment", df_name = "sigUp_Adama_U"),
      list(analysis = "GO_enrichment", df_name = "sigUp_Erer_U"),
      list(analysis = "GO_enrichment", df_name = "sigUp_Jijiga_U")
    ),
    title = "GO Enrichment Terms shared among sites",
    filename = "Venn_Sites_GO_enrichment.png"
  ),
  
  # --- A 3-set Venn for Up-regulated Wild Caught Larvae GSEA ---
  sites_GSEA = list(
    dfs_to_compare = list(
      list(analysis = "GSEA", df_name = "Adama_vs_Erer_GSEA_results.csv"),
      list(analysis = "GSEA", df_name = "Adama_vs_Jijiga_GSEA_results.csv"),
      list(analysis = "GSEA", df_name = "Erer_vs_Jijiga_GSEA_results.csv")
    ),
    title = "GSEA Terms shared among sites",
    filename = "Venn_Sites_GSEA.png"
  )
  
)


# --- 2. Loop Through Definitions and Generate Diagrams ---

# This loop iterates through each set of instructions you defined above.
for (venn_name in names(venn_definitions)) {
  
  current_venn <- venn_definitions[[venn_name]]
  
  # Create an empty list to hold the vectors of IDs for the diagram
  list_of_id_sets <- list()
  
  # Create a vector to hold the names for each circle in the diagram
  set_names <- c()
  
  # Extract the IDs from each specified data frame
  for (df_info in current_venn$dfs_to_compare) {
    
    # Safely get the data frame from the nested list structure
    df <- all_data_nested[[df_info$analysis]][[df_info$df_name]]
    
    if (!is.null(df)) {
      # Add the vector of unique IDs to our list
      list_of_id_sets <- append(list_of_id_sets, list(unique(df$ID)))
      
      # Add the name for this set (stripping .csv for cleaner labels)
      set_names <- c(set_names, sub("\\.csv$", "", df_info$df_name))
      
    } else {
      warning(paste("Data frame not found:", df_info$df_name, "in analysis", df_info$analysis, ". Skipping for this Venn."))
    }
  }
  
  # Check if we have a valid number of sets to compare
  if (length(list_of_id_sets) < 2 || length(list_of_id_sets) > 5) {
    warning(paste("Skipping '", venn_name, "' because it has ", length(list_of_id_sets), " sets. The package supports 2-5 sets.", sep=""))
    next # Skip to the next iteration of the loop
  }
  
  # Name the list of ID sets. This is required by the venn.diagram function.
  names(list_of_id_sets) <- set_names
  
  # --- Generate the Venn Diagram and Save to File ---
  
  # Suppress the default log file that the package creates
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  
  venn.diagram(
    x = list_of_id_sets,
    filename = current_venn$filename,
    
    # --- Output settings ---
    imagetype = "png", # You can also use "tiff"
    height = 1000,
    width = 1200,
    resolution = 600,
    compression = "lzw",
    
    # --- Title and Labels ---
    main = current_venn$title,
    main.cex = 0.6,
    sub.cex = 0.5,
    
    # --- Circles ---
    # Set the colors for the circles
    fill = c("dodgerblue", "goldenrod1", "darkorange1", "violetred1", "seagreen3")[1:length(list_of_id_sets)],
    col = "transparent",
    alpha = 0.50,
    
    # --- Numbers inside circles ---
    cex = 0.5,
    fontface = "bold",
    
    # --- Set Names (Labels for each circle) ---
    category.names = gsub("sigUp_|sigDown_|_U|_D|_GSEA_results|\\.csv", "", set_names),
    cat.cex = 0.4,
    cat.fontface = "bold",
    cat.dist = 0.05,
    cat.pos = 0
  )
  
  cat(paste("Successfully generated Venn diagram and saved to '", current_venn$filename, "'\n", sep=""))
}

# Clean up variables
rm(venn_definitions, current_venn, list_of_id_sets, set_names, df_info, df, venn_name)







