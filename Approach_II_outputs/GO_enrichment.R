# R Script for Non-Directional Filtering and GO Enrichment Analysis
# This version runs the analysis on BOTH the unfiltered and filtered gene sets.

# Setting working directory --------

setwd("C:/Users/aaisu/Box/Carter Lab/Projects/RNAseq_EthiopiaNUCI/Approach_II_outputs")

# Consolidated GO Enrichment Pipeline for Unfiltered and Directionally Filtered Gene Sets

# --- 0. LOAD LIBRARIES and ASSUMPTIONS ---

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
# Ensure your custom annotation DB packages are loaded
library(org.Astephensi.eg.db)
library(GO.db)

# This script assumes you have:
#   1. A folder containing the full DESeq2 result CSVs.
#   2. Your custom annotation objects (ann, term_names) or the code to generate them.

# --- 1. DEFINE INPUT AND OUTPUT DIRECTORIES ---

results_directory <- "./DESeq2_full_results/"

# Define a clear, hierarchical folder structure for all outputs
dir.create("./GOenrichment/unfiltered_results/", showWarnings = FALSE, recursive = TRUE)
dir.create("./GOenrichment/unfiltered_plots/", showWarnings = FALSE, recursive = TRUE)
dir.create("./GOenrichment/directional_filtered_results/", showWarnings = FALSE, recursive = TRUE)
dir.create("./GOenrichment/directional_filtered_plots/", showWarnings = FALSE, recursive = TRUE)


# --- 2. CREATE DIRECTIONAL MASTER LISTS of STAGE GENES and GENE UNIVERSE ---

cat("--- Step 1: Identifying UP/DOWN stage genes and defining the universe ---\n")
all_results_files <- list.files(path = results_directory, pattern = "\\.csv$", full.names = TRUE)

stage_genes_UP <- c()
stage_genes_DOWN <- c()
universe_genes <- c()

for (csv_file in all_results_files) {
  results_df <- read.csv(csv_file)
  universe_genes <- c(universe_genes, results_df$geneID)
  
  if (grepl("IS2_vs_IS1|IS3_vs_IS2|IS4_vs_IS3", basename(csv_file))) {
    cat("Reading stage-specific genes from:", basename(csv_file), "\n")
    # Get UP-regulated stage genes
    sig_up <- subset(results_df, padj < 0.05 & log2FoldChange > 1)$geneID
    stage_genes_UP <- c(stage_genes_UP, sig_up)
    
    # Get DOWN-regulated stage genes
    sig_down <- subset(results_df, padj < 0.05 & log2FoldChange < -1)$geneID
    stage_genes_DOWN <- c(stage_genes_DOWN, sig_down)
  }
}

# Keep only unique gene IDs in each list
stage_genes_UP <- unique(stage_genes_UP)
stage_genes_DOWN <- unique(stage_genes_DOWN)
universe_genes <- unique(universe_genes)

cat("... Identified", length(stage_genes_UP), "unique UP-regulated stage genes for filtering.\n")
cat("... Identified", length(stage_genes_DOWN), "unique DOWN-regulated stage genes for filtering.\n")
cat("... Defined a universe of", length(universe_genes), "total tested genes.\n")


# --- (Optional) Annotation DB Preparation ---
# It's good practice to have this here if not loaded from an RDS file
# Before running GSEA taking a look at keytypes in org.db
keytypes(org.Astephensi.eg.db)

# Set the org.db
orgdb <- org.Astephensi.eg.db

# Pull GO annotations (direct only) and keep BP
ann <- AnnotationDbi::select(
  orgdb,
  keys = keys(orgdb, keytype = "GID"), # or whatever ID you want to use
  columns = c("GO","ONTOLOGY"),
  keytype = "GID"
)
ann <- ann[!is.na(ann$GO), c("GO","GID")] # Removed filter for only Biological processes -  & ann$ONTOLOGY == "BP"
colnames(ann) <- c("term", "gene")  # TERM2GENE

# 2) Optional TERM2NAME (GO term names)
term_names <- AnnotationDbi::select(GO.db, keys = unique(ann$term),
                                    columns = "TERM", keytype = "GOID")
colnames(term_names) <- c("term","name")  # TERM2NAME



# --- 3. CONSOLIDATED GO ENRICHMENT LOOP ---

cat("\n--- Step 2: Running GO Enrichment on all gene sets ---\n")

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
      res_dir <- "./GOenrichment/unfiltered_results/"
      plt_dir <- "./GOenrichment/unfiltered_plots/"
      suffix <- "unfiltered"
    } else {
      res_dir <- "./GOenrichment/directional_filtered_results/"
      plt_dir <- "./GOenrichment/directional_filtered_plots/"
      suffix <- "dir_filtered"
    }
    
    # Run GO enrichment with error handling
    tryCatch({
      enricher_res <- enricher(
        gene = gene_list, TERM2GENE = ann, TERM2NAME = term_names,
        universe = universe_genes, pAdjustMethod = "BH", pvalueCutoff = 0.05,
        qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500
      )
      
      if (is.null(enricher_res) || nrow(enricher_res@result) == 0) {
        cat("  -> No significant GO terms found.\n")
      } else {
        # Save results
        file_key <- paste0(ifelse(current_analysis$direction == "Up", "sigUp", "sigDown"), "_", base_name)
        write.csv(as.data.frame(enricher_res@result),
                  file = file.path(res_dir, paste0(file_key, "_GO_", suffix, ".csv")),
                  row.names = FALSE)
        
        # Save plots
        plot_title <- paste("Top GO Terms (", suffix, ") for\n", file_key)
        p1 <- barplot(enricher_res, showCategory = 20, title = plot_title)
        ggsave(filename = file.path(plt_dir, paste0(file_key, "_GO_barplot_", suffix, ".png")),
               plot = p1, width = 10, height = 10)
        
        p2 <- dotplot(enricher_res, showCategory = 20, title = plot_title)
        ggsave(filename = file.path(plt_dir, paste0(file_key, "_GO_dotplot_", suffix, ".png")),
               plot = p2, width = 10, height = 10)
        
        cat("  -> Successfully saved", suffix, "GO results and plots.\n")
      }
    }, error = function(e) {
      cat("  -> ERROR during GO enrichment analysis.\n")
      print(e)
    })
  }
}

cat("\n--- Consolidated GO Enrichment Analysis Pipeline Complete ---\n")