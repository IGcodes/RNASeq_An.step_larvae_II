# R Script to run KEGG Pathway Analysis on significant genes from saved DESeq2 results.

# Setting working directory
setwd("C:/Users/aaisu/Box/Carter Lab/Projects/RNAseq_EthiopiaNUCI/Approach_II_outputs")


# Consolidated KEGG Analysis Pipeline for Unfiltered and Directionally Filtered Gene Sets

# --- 0. LOAD LIBRARIES and ASSUMPTIONS ---

# This script assumes you have:
#   1. The 'clusterProfiler' and 'enrichplot' packages installed.
#   2. A folder containing the full DESeq2 result CSVs.
#   3. Your gene IDs are in a format compatible with enrichKEGG (e.g., NCBI Gene IDs).

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# --- 1. DEFINE INPUT AND OUTPUT DIRECTORIES ---

results_directory <- "./DESeq2_full_results/"

# Define a clear, hierarchical folder structure for all outputs
dir.create("./KEGGS/unfiltered_results/", showWarnings = FALSE, recursive = TRUE)
dir.create("./KEGGS/unfiltered_plots/", showWarnings = FALSE, recursive = TRUE)
dir.create("./KEGGS/directional_filtered_results/", showWarnings = FALSE, recursive = TRUE)
dir.create("./KEGGS/directional_filtered_plots/", showWarnings = FALSE, recursive = TRUE)


# --- 2. CREATE DIRECTIONAL MASTER LISTS OF STAGE-RELATED GENES ---

cat("--- Step 1: Identifying UP and DOWN stage-related genes for filtering ---\n")
all_results_files <- list.files(path = results_directory, pattern = "\\.csv$", full.names = TRUE)

stage_genes_UP <- c()
stage_genes_DOWN <- c()

for (csv_file in all_results_files) {
  if (grepl("IS2_vs_IS1|IS3_vs_IS2|IS4_vs_IS3", basename(csv_file))) {
    cat("Reading stage-specific genes from:", basename(csv_file), "\n")
    stage_results_df <- read.csv(csv_file)
    
    # Get significant UP genes and add to the UP list
    sig_up <- subset(stage_results_df, padj < 0.05 & log2FoldChange > 1)$geneID
    stage_genes_UP <- c(stage_genes_UP, sig_up)
    
    # Get significant DOWN genes and add to the DOWN list
    sig_down <- subset(stage_results_df, padj < 0.05 & log2FoldChange < -1)$geneID
    stage_genes_DOWN <- c(stage_genes_DOWN, sig_down)
  }
}

# Keep only unique gene IDs in each list
stage_genes_UP <- unique(stage_genes_UP)
stage_genes_DOWN <- unique(stage_genes_DOWN)
cat("... Identified", length(stage_genes_UP), "unique UP-regulated stage genes.\n")
cat("... Identified", length(stage_genes_DOWN), "unique DOWN-regulated stage genes.\n")


# --- 3. CONSOLIDATED KEGG ANALYSIS LOOP ---

cat("\n--- Step 2: Running KEGG analysis on all gene sets ---\n")

processed_gene_lists <- c()

for (csv_file in all_results_files) {
  # Skip the stage comparison files themselves for the main analysis
  # if (grepl("IS2_vs_IS1|IS3_vs_IS2|IS4_vs_IS3", basename(csv_file))) { next }
  
  # --- Data Preparation for the current comparison ---
  base_name <- sub("\\.csv$", "", sub("results_", "", basename(csv_file)))
  groups <- strsplit(base_name, "_vs_")[[1]]
  full_results_df <- read.csv(csv_file)
  
  # Get unfiltered lists
  sigUp_unfiltered <- subset(full_results_df, padj < 0.05 & log2FoldChange > 1)$geneID
  sigDown_unfiltered <- subset(full_results_df, padj < 0.05 & log2FoldChange < -1)$geneID
  
  # Get directionally filtered lists
  sigUp_filtered <- setdiff(sigUp_unfiltered, stage_genes_UP)
  sigDown_filtered <- setdiff(sigDown_unfiltered, stage_genes_DOWN)
  
  # --- Define all analyses to run for this one comparison file ---
  analyses <- list(
    unfiltered_sigUp = list(genes = sigUp_unfiltered, type = "unfiltered", direction = "Up"),
    unfiltered_sigDown = list(genes = sigDown_unfiltered, type = "unfiltered", direction = "Down"),
    filtered_sigUp = list(genes = sigUp_filtered, type = "filtered", direction = "Up"),
    filtered_sigDown = list(genes = sigDown_filtered, type = "filtered", direction = "Down")
  )
  
  # --- Inner loop to execute the four defined analyses ---
  for (analysis_name in names(analyses)) {
    
    current_analysis <- analyses[[analysis_name]]
    gene_list <- current_analysis$genes
    
    # Check for redundancy
    reciprocal_key <- paste(
      current_analysis$type, 
      ifelse(current_analysis$direction == "Up", "Down", "Up"), 
      groups[2], groups[1], sep = "_"
    )
    current_key <- paste(current_analysis$type, current_analysis$direction, groups[1], groups[2], sep = "_")
    
    if (reciprocal_key %in% processed_gene_lists) { next }
    processed_gene_lists <- c(processed_gene_lists, current_key)
    
    cat("\n--------------------------------------------------\n")
    cat("Processing:", current_key, "\n")
    cat("  -> Number of genes:", length(gene_list), "\n")

    if (length(gene_list) == 0) {
      cat("  -> No genes in this list. Skipping analysis.\n")
      next
    }
    
    # Set output paths and suffixes based on analysis type
    if (current_analysis$type == "unfiltered") {
      res_dir <- "./KEGGS/unfiltered_results/"
      plt_dir <- "./KEGGS/unfiltered_plots/"
      suffix <- "unfiltered"
    } else {
      res_dir <- "./KEGGS/directional_filtered_results/"
      plt_dir <- "./KEGGS/directional_filtered_plots/"
      suffix <- "dir_filtered"
    }
    
    # Run KEGG analysis with error handling
    tryCatch({
      kegg_res <- enrichKEGG(
        gene = gene_list, organism = "aste", keyType = "ncbi-geneid",
        pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.10
      )
      
      if (is.null(kegg_res) || nrow(kegg_res@result) == 0) {
        cat("  -> No significant KEGG pathways found.\n")
      } else {
        # Save results
        file_key <- paste0(current_analysis$direction, "_", base_name)
        write.csv(as.data.frame(kegg_res@result),
                  file = file.path(res_dir, paste0(file_key, "_KEGG_", suffix, ".csv")),
                  row.names = FALSE)
        
        # Save plots
        plot_title <- paste("Top KEGG Pathways (", suffix, ") for\n", file_key)
        p1 <- barplot(kegg_res, showCategory = 20, title = plot_title)
        ggsave(filename = file.path(plt_dir, paste0(file_key, "_KEGG_barplot_", suffix, ".png")),
               plot = p1, width = 10, height = 10)
        
        p2 <- dotplot(kegg_res, showCategory = 20, title = plot_title)
        ggsave(filename = file.path(plt_dir, paste0(file_key, "_KEGG_dotplot_", suffix, ".png")),
               plot = p2, width = 10, height = 10)
               
        cat("  -> Successfully saved", suffix, "KEGG results and plots.\n")
      }
    }, error = function(e) {
      cat("  -> ERROR during enrichKEGG analysis. Likely no genes mapped to KEGG.\n")
      print(e)
    })
  }
}

cat("\n--- Consolidated KEGG Analysis Pipeline Complete ---\n")