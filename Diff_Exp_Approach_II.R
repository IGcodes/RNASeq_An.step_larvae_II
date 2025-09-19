# Setting path to the working directory
setwd("C:/Users/aaisu/Box/Carter Lab/Projects/RNAseq_EthiopiaNUCI")
## Written on R version 4.5.1

# Importing necessary libraries
library(rtracklayer)
library(GenomicFeatures)
library(DESeq2)
library(dplyr)
library(tibble)
library(ggplot2)
library(tximport)
library(readr)
library(AnnotationDbi)     # if you need to map IDs
library(ashr)            # for LFC shrinkage
library(pheatmap)
library(clusterProfiler)
library(enrichplot)
library(AnnotationForge)
library(txdbmaker)
library(biomaRt)
library(org.Astephensi.eg.db)
library(GO.db)
library(ggrepel)
library(patchwork)
library(RColorBrewer)

#  --- 1. Read in sample metadata
# Expect a table with columns: sampleName, stage, quant_dir
sampleTable <- read_csv("./quant_files/All_sample_info.csv")
# e.g. sample_info.csv:
# sampleName,stage,quant_dir
# S1,L1,/path/to/S1
# S2,L1,/path/to/S2
# ...
sampleTable$stage <- factor(sampleTable$stage, levels=c("IS1","IS2","IS3","IS4", "DK"))
sampleTable$Collection <- factor(sampleTable$Collection, levels = c("Lab", "Wild"))
sampleTable$Site <- factor(sampleTable$Site, levels = c("UCI", "Adama", "Erer", "Jijiga") )
rownames(sampleTable) <- sampleTable$sampleName

#  --- 2. Build vector of Salmon quant.sf files
Collection_files <- file.path(sampleTable$quant_dir, "quant.sf")
names(Collection_files) <- sampleTable$sampleName

Lab_files <- file.path(sampleTable[which(sampleTable$Collection == "Lab"),]$quant_dir, "quant.sf")
names(Lab_files) <- sampleTable[which(sampleTable$Collection == "Lab"),]$sampleName

Wild_files <- file.path(sampleTable[which(sampleTable$Collection == "Wild"),]$quant_dir, "quant.sf")
names(Wild_files) <- sampleTable[which(sampleTable$Collection == "Wild"),]$sampleName

#  --- 3. Import transcript-level counts and aggregate to genes
# You need a two-column data.frame tx2gene mapping transcripts to genes
transcript_info <- read.csv("./Reference_transcriptome/Anstep_UCI_V1.0_transcripts_info.csv")
# Each transcript ID ends with ".1" extension indicating the version. However, this has been disregarded in the Salmon quantification.
# Therefore I had to remove it from the tx2gene table
# Note: Remember to add that part before any conversion of transcript IDs to another format.
# For the KEGGS analysis we need the gene ID without the perfix "LOC". Therefore it was also removed
transcript_info$accession <- sub("\\.1$", "", transcript_info$accession)
transcript_info$gene_id <- sub("^LOC", "", transcript_info$gene_id)
tx2gene <- transcript_info[,c(1,5)]

# Reading in files
Collection_txi <- tximport(Collection_files, type="salmon", tx2gene=tx2gene, ignoreTxVersion=TRUE)
Lab_txi <- tximport(Lab_files, type="salmon", tx2gene=tx2gene, ignoreTxVersion=TRUE)
Wild_txi <- tximport(Wild_files, type="salmon", tx2gene=tx2gene, ignoreTxVersion=TRUE)

# --- 4. CREATE CONSOLIDATED DESEQ DATASETS ---
# Instead of many dds objects, we only need one for each main experimental design.
# All comparisons for a given factor will be extracted from the SAME object.

dds_stage <- DESeqDataSetFromTximport(Lab_txi, colData = sampleTable[which(sampleTable$Collection == "Lab"),], design = ~ stage)
dds_site <- DESeqDataSetFromTximport(Wild_txi, colData = sampleTable[which(sampleTable$Collection == "Wild"),], design = ~ Site)
dds_collection <- DESeqDataSetFromTximport(Collection_txi, colData = sampleTable, design = ~ Collection)


# --- 5. PREFILTER ONCE ---
# We create the 'keep' vector once and apply it to our consolidated objects.
smallestGroupSize <- 3
keep_stage <- rowSums(counts(dds_stage) >= 10) >= smallestGroupSize
keep_site <- rowSums(counts(dds_site) >= 10) >= smallestGroupSize
keep_collection <- rowSums(counts(dds_collection) >= 10) >= smallestGroupSize

dds_stage <- dds_stage[keep_stage,]
dds_site <- dds_site[keep_site,]
dds_collection <- dds_collection[keep_collection,]


# --- 6. RUN THE DESEQ PIPELINE ONCE PER DESIGN ---
# This is the main computational step. It only needs to be run once for each design.
cat("Running DESeq on the 'stage' model...\n")
dds_stage <- DESeq(dds_stage)

cat("Running DESeq on the 'Site' model...\n")
dds_site <- DESeq(dds_site)

cat("Running DESeq on the 'Collection' model...\n")
dds_collection <- DESeq(dds_collection)


# --- 7. DEFINE ALL COMPARISONS IN A "PLAN" ---
# This is the core of the automation. We create a data frame where each row
# is a comparison we want to perform. This is much easier to manage and extend.

# Create a dedicated folder to store the results CSVs.
output_directory <- "./Approach_II_outputs/DESeq2_full_results/"
dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)


comparisons_plan <- data.frame(
  # The DESeqDataSet object to use
  dds_object_name = c("dds_stage", "dds_stage", "dds_stage", "dds_site", "dds_site", "dds_site", "dds_site", "dds_site", "dds_site", "dds_collection"),
  # The factor to test
  factor = c("stage", "stage", "stage", "Site", "Site", "Site", "Site", "Site", "Site", "Collection"),
  # Group 1 (numerator)
  group1 = c("IS2", "IS3", "IS4", "Erer", "Adama", "Jijiga", "Adama", "Jijiga", "Erer", "Wild"),
  # Group 2 (denominator/reference)
  group2 = c("IS1", "IS2", "IS3", "Jijiga", "Jijiga", "Erer", "Erer", "Adama", "Adama", "Lab")
)


# --- 8. AUTOMATED LOOP FOR RESULTS EXTRACTION ---
# This single loop replaces all of your repetitive lfcShrink and subset calls.

# Create an empty list to store all our results data frames
all_results <- list()

for (i in 1:nrow(comparisons_plan)) {
  
  # Get details for the current comparison from our plan
  plan_row <- comparisons_plan[i, ]
  dds_obj <- get(plan_row$dds_object_name) # 'get()' retrieves the R object by its string name
  
  # Define the contrast vector for lfcShrink
  current_contrast <- c(plan_row$factor, plan_row$group1, plan_row$group2)
  
  cat("Processing:", paste(current_contrast, collapse=" "), "\n")
  
  # Run lfcShrink using the contrast argument
  resLFC <- lfcShrink(dds_obj, contrast = current_contrast, type = "ashr")
  
  # --- NEW: SAVE FULL RESULTS TO CSV ---
  
  # 1. Convert the complete results object to a data frame
  full_results_df <- as.data.frame(resLFC)
  
  # 2. Add the gene IDs (which are stored as row names) as a new column
  full_results_df$geneID <- rownames(full_results_df)
  
  # 3. Reorder columns to make geneID the first column for readability
  full_results_df <- full_results_df[, c("geneID", setdiff(names(full_results_df), "geneID"))]
  
  # 4. Create a descriptive filename
  output_filename <- paste0("results_", plan_row$group1, "_vs_", plan_row$group2, ".csv")
  full_output_path <- file.path(output_directory, output_filename)
  
  # 5. Write the data frame to the CSV file
  write.csv(full_results_df, file = full_output_path, row.names = FALSE)
  cat("  -> Full results saved to:", full_output_path, "\n")
  
  # --- END OF NEW SECTION ---
  
  # Create a descriptive name for the results
  res_name <- paste0("resLFC_", plan_row$group1, "_vs_", plan_row$group2)
  
  # Subset for significantly up- and down-regulated genes
  sigUp <- subset(resLFC, padj < 0.05 & log2FoldChange > 1)
  sigDown <- subset(resLFC, padj < 0.05 & log2FoldChange < -1)
  
  # Store the results in our list with clean names
  all_results[[paste0("sigUp_", plan_row$group1, "_vs_", plan_row$group2)]] <- as.data.frame(sigUp)
  all_results[[paste0("sigDown_", plan_row$group1, "_vs_", plan_row$group2)]] <- as.data.frame(sigDown)
}

cat("\n--- Analysis Complete ---\n")


# --- 9. ACCESSING RESULTS ---
# All your results are now neatly organized in a single list.
# You can see all the result sets you created:
# names(all_results)

# You can access a specific result data frame like this:
# head(all_results$sigUp_IS2_vs_IS1)

#### Testing reciprocity of the gene lists as way of sanity check ################
# --- 1. IDENTIFY ALL UP-REGULATED RESULT SETS ---

# Get all the names from your results list
all_result_names <- names(all_results)

# Filter this list to get only the names of the "sigUp" data frames
up_regulated_names <- all_result_names[startsWith(all_result_names, "sigUp_")]

cat("Found", length(up_regulated_names), "up-regulated result sets to test.\n\n")


# --- 2. CREATE A DATA FRAME TO STORE TEST RESULTS ---

# This data frame will give us a neat summary of the validation checks.
reciprocity_summary <- data.frame(
  UpComparison = character(),
  DownReciprocal = character(),
  AreGeneListsIdentical = logical(),
  NumberOfGenes = integer(),
  stringsAsFactors = FALSE
)


# --- 3. LOOP THROUGH AND TEST EACH PAIR ---

for (up_name in up_regulated_names) {
  
  # 3a. Programmatically determine the name of the reciprocal down-regulated set
  # Example: up_name = "sigUp_Jijiga_vs_Erer"
  
  # Remove the "sigUp_" prefix to get the base comparison
  base_comparison <- sub("sigUp_", "", up_name) # -> "Jijiga_vs_Erer"
  
  # Split the base comparison into its two parts
  parts <- strsplit(base_comparison, "_vs_")[[1]] # -> c("Jijiga", "Erer")
  
  # Reverse the parts and paste them back together
  reciprocal_base <- paste(parts[2], parts[1], sep = "_vs_") # -> "Erer_vs_Jijiga"
  
  # Prepend "sigDown_" to get the final name of the reciprocal data frame
  down_name <- paste0("sigDown_", reciprocal_base) # -> "sigDown_Erer_vs_Jijiga"
  
  # 3b. Check if the expected reciprocal result actually exists
  if (!down_name %in% all_result_names) {
    cat("Warning: Could not find the reciprocal partner for", up_name, ". Skipping.\n")
    next # Skip to the next iteration
  }
  
  # 3c. Retrieve the two data frames from the main list
  up_regulated_df <- all_results[[up_name]]
  down_regulated_df <- all_results[[down_name]]
  
  # 3d. Extract the gene IDs.
  # DESeq2 results store gene IDs as row names.
  up_regulated_genes <- rownames(up_regulated_df)
  down_regulated_genes <- rownames(down_regulated_df)
  
  # 3e. Perform the test
  # The identical() function checks if two R objects are exactly equal,
  # including the order of elements, which is what we want here.
  are_identical <- identical(up_regulated_genes, down_regulated_genes)
  
  # 3f. Add the results to our summary data frame
  reciprocity_summary <- rbind(reciprocity_summary, data.frame(
    UpComparison = up_name,
    DownReciprocal = down_name,
    AreGeneListsIdentical = are_identical,
    NumberOfGenes = length(up_regulated_genes)
  ))
}


# --- 4. DISPLAY THE FINAL SUMMARY ---

cat("\n--- Reciprocity Test Summary ---\n")
print(reciprocity_summary)


# --- 10. Plotting PCA wiht VST ---

# Creating vst objects
vsd_stage <- vst(dds_stage, blind=FALSE)
vsd_site <- vst(dds_site, blind=FALSE)
vsd_collection <- vst(dds_collection, blind=FALSE)


# --- 2. DEFINE CONSISTENT AESTHETICS ---
site_colors <- brewer.pal(4, "Set2")
names(site_colors) <- c("UCI", "Adama", "Erer", "Jijiga")

stage_colors <- brewer.pal(4, "Dark2")
names(stage_colors) <- c("IS1", "IS2", "IS3", "IS4")

stage_shapes <- c("IS1" = 16, "IS2" = 17, "IS3" = 15, "IS4" = 18)


# --- 3. REBUILD PLOT A (PCA of Lab Samples by Stage) ---

pca_data_stage <- plotPCA(vsd_stage, intgroup = "stage", returnData = TRUE)
percentVar_stage <- round(100 * attr(pca_data_stage, "percentVar"))

stage_PCA_plot <- ggplot(pca_data_stage, aes(x = PC1, y = PC2, color = stage, shape = stage)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(aes(label = name), size = 4.5, max.overlaps = Inf, box.padding = 0.5) +
  scale_color_manual(values = stage_colors) +
  scale_shape_manual(values = stage_shapes) +
  xlab(paste0("PC1: ", percentVar_stage[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_stage[2], "% variance")) +
  ggtitle("PCA of Lab Samples by Stage") +
  theme_bw(base_size = 16)


# --- 4. REBUILD PLOT B (PCA of Wild Samples by Site) ---

pca_data_site <- plotPCA(vsd_site, intgroup = "Site", returnData = TRUE)
percentVar_site <- round(100 * attr(pca_data_site, "percentVar"))

site_PCA_plot <- ggplot(pca_data_site, aes(x = PC1, y = PC2, color = Site)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(aes(label = name), size = 4.5, max.overlaps = Inf, box.padding = 0.5) +
  stat_ellipse(type = "norm", level = 0.95, linewidth = 1) +
  scale_color_manual(values = site_colors) +
  xlab(paste0("PC1: ", percentVar_site[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_site[2], "% variance")) +
  ggtitle("PCA of Wild Samples by Site") +
  theme_bw(base_size = 16)


# --- 5. REBUILD PLOT C (PCA of All Samples) ---

pca_data_collection <- plotPCA(vsd_collection, intgroup = c("Collection", "Site"), returnData = TRUE)
percentVar_collection <- round(100 * attr(pca_data_collection, "percentVar"))

collection_PCA_plot <- ggplot(pca_data_collection, aes(x = PC1, y = PC2, color = Site, shape = Collection)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = Site), type = "norm", level = 0.95, linewidth = 1) +
  scale_color_manual(values = site_colors) +
  scale_shape_manual(values = c("Lab" = 16, "Wild" = 17)) +
  xlab(paste0("PC1: ", percentVar_collection[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_collection[2], "% variance")) +
  ggtitle("PCA of All Samples") +
  theme_bw(base_size = 16)


# --- 6. COMBINE, ANNOTATE, AND SAVE THE FINAL FIGURE ---

# 6a. *** LAYOUT CHANGE: Stack all three plots vertically ***
# The '/' operator in patchwork stacks plots.
final_layout <- stage_PCA_plot + site_PCA_plot + collection_PCA_plot

# 6b. Assemble the final plot with annotations
final_plot <- final_layout +
  plot_layout(guides = "collect") + # Collect legends into one area
  plot_annotation(
    title = 'Transcriptomic Variation in An. stephensi Larvae',
    tag_levels = 'a',
    theme = theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5))
  ) & theme(legend.position = "bottom", legend.box = "vertical")


# 6c. Save the final plot with dimensions adjusted for a vertical layout
output_filename <- "./Approach_II_outputs/PCA_plots/Anstep_Larvae_PCA_plots_horizontal.png"
dir.create(dirname(output_filename), showWarnings = FALSE, recursive = TRUE)

# Adjusted width and height for a tall, single-column figure
ggsave(output_filename, plot = final_plot, width = 16, height = 8, units = "in", dpi = 300)

cat("Final combined plot saved to:", output_filename, "\n")



####################################################################################################################################
################# GSEA for all comparisons ####################################################################################

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


# --- 1. SETUP OUTPUT DIRECTORIES ---

# We create the directories where the results will be saved.
dir.create("./Approach_II_outputs/GSEA/All_samples/Dataframes", showWarnings = FALSE)
dir.create("./Approach_II_outputs/GSEA/All_samples/Graphs/", showWarnings = FALSE)


# --- 2. AUTOMATED GSEA LOOP ---

cat("Starting automated GSEA analysis...\n")

# This vector will keep track of the comparisons we've already processed
# to avoid running the reciprocal (e.g., B vs A after A vs B).
processed_pairs <- c()

# Loop through each comparison defined in your plan
for (i in 1:nrow(comparisons_plan)) {
  
  # 2a. Get details for the current comparison
  plan_row <- comparisons_plan[i, ]
  
  # 2b. Check for and skip reciprocal pairs
  # We create a unique, sorted "canonical" name for each pair.
  # For both "Erer vs. Jijiga" and "Jijiga vs. Erer", this will be "Erer_Jijiga".
  pair_name <- paste(sort(c(plan_row$group1, plan_row$group2)), collapse = "_")
  
  if (pair_name %in% processed_pairs) {
    cat("\nSkipping reciprocal comparison:", plan_row$group1, "vs.", plan_row$group2, "\n")
    next # Skip to the next iteration
  }
  
  # If it's a new pair, add it to our list of processed pairs and proceed.
  processed_pairs <- c(processed_pairs, pair_name)
  cat("\n--------------------------------------------------\n")
  cat("Processing:", plan_row$group1, "vs.", plan_row$group2, "\n")
  
  # 2c. Get the full DESeq2 results table for this comparison
  # This provides the complete, unshrunk results needed for ranking.
  dds_obj <- get(plan_row$dds_object_name)
  
  # Note: GSEA works best on unshrunk LFCs. We use results() here.
  # If you still prefer shrunken values, you can use lfcShrink() as before.
  full_results <- results(dds_obj, contrast = c(plan_row$factor, plan_row$group1, plan_row$group2))
  
  # 2d. Prepare the ranked gene list for GSEA
  # GSEA requires a named vector of numeric values, sorted in descending order.
  
  gene_list <- full_results$log2FoldChange
  names(gene_list) <- rownames(full_results)
  
  # Remove any genes that have NA values, which can cause errors
  gene_list <- gene_list[!is.na(gene_list)]
  
  # Sort the list in descending order (most up-regulated to most down-regulated)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # 2e. Run GSEA
  cat("  -> Running GSEA...\n")
  gsea_result <- GSEA(
    geneList = gene_list,
    TERM2GENE = ann,
    TERM2NAME = term_names,
    pvalueCutoff = 0.05,
    minGSSize = 10,
    maxGSSize = 500,
    verbose = FALSE # Set to TRUE if you want more detailed progress output
  )
  
  # 2f. Save results and plots
  # Create a clean and descriptive file identifier
  file_identifier <- paste0(plan_row$group1, "_vs_", plan_row$group2)
  
  # Check if GSEA returned any significant results
  if (is.null(gsea_result) || nrow(gsea_result@result) == 0) {
    cat("  -> No significant GSEA terms found for this comparison.\n")
  } else {
    # Save the GSEA results table to a CSV file
    gsea_results_df <- as.data.frame(gsea_result@result)
    write.csv(gsea_results_df, file = file.path("./Approach_II_outputs/GSEA/All_samples/Dataframes", paste0(file_identifier, "_GSEA_results.csv")), row.names = FALSE)
    
    # Generate and save the Dotplot
    p1 <- dotplot(gsea_result, showCategory = 20, title = paste(file_identifier, "GSEA Dotplot"))
    ggsave(filename = file.path("./Approach_II_outputs/GSEA/All_samples/Graphs/", paste0(file_identifier, "_GSEA_dotplot.png")), plot = p1, width = 10, height = 10)
    
    # Generate and save the Ridgeplot
    p2 <- ridgeplot(gsea_result) + labs(title = paste(file_identifier, "GSEA Ridgeplot"))
    ggsave(filename = file.path("./Approach_II_outputs/GSEA/All_samples/Graphs/", paste0(file_identifier, "_GSEA_ridgeplot.png")), plot = p2, width = 10, height = 12)
    
    cat("  -> Successfully saved results and plots.\n")
  }
}

cat("\n--- GSEA Pipeline Complete ---\n")
