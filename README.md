# Transcriptomic Analysis of <i> Anopheles stephensi </i> Larvae

This repository contains a complete bioinformatics pipeline for the analysis of RNA-seq data from <i> Anopheles stephensi </i> larvae. The workflow covers every step from reference transcriptome preparation and read quantification to differential gene expression, functional enrichment analysis (GO, KEGG, GSEA), and downstream comparative analysis of pathway results.

The primary goal of this analysis is to compare gene expression profiles across different larval stages (IS1, IS2, IS3, IS4) and between wild-caught and lab-reared mosquito populations from various sites (UCI, Adama, Erer, Jijiga).

## 🧬 Analysis Workflow
The analysis is designed to be run sequentially, where the outputs of earlier steps become the inputs for later ones.

1. Reference & Annotation Preparation (```/Reference_transcriptome```)

    * Generate a transcript-to-gene mapping file from the reference FASTA using ```tx2gene_creator.py```.
  
    * Build a custom R ```OrgDb``` annotation package from a Gene Ontology Annotation (GAF) file using ```gaf2orgDB.R```.

2. Read Quantification (```/shell_scripts```)

    * Index the reference transcriptome for quasi-mapping using ```indexing_reference_salmon.sh```.
  
    * Quantify transcript abundance from paired-end FASTQ files for each sample using ```quant_salmon.sh```.

3. Core Differential Expression & GSEA (```/Diff_Exp_Approach_II.R```)

    * Import Salmon quantification data using ```tximport```.
    
    * Perform differential expression analysis between specified sample groups using ```DESeq2```.
    
    * Generate PCA plots to visualize sample clustering.
    
    * Run Gene Set Enrichment Analysis (GSEA) on all comparisons.

4. Functional Enrichment Analysis (```/Approach_II_outputs```)

    * Perform KEGG pathway over-representation analysis on significant gene lists using KEGGS_analysis.R.
    
    * Perform GO term over-representation analysis on significant gene lists using GO_enrichment.R.

5. Downstream Pathway Comparison (```/Approach_II_outputs/Pathway_selection```)

    * Integrate, filter, and compare the results from the various enrichment analyses to identify shared and unique pathways using Pathway_selections.R.


## 📂 Repository Structure

.

├── Diff_Exp_Approach_II.R

├── Approach_II_outputs

│   ├── GO_enrichment.R

│   ├── KEGGS_analysis.R

│   └── Pathway_selection

│       └── Pathway_selections.R

├── Reference_transcriptome

│   ├── gaf2orgDB.R

│   └── tx2gene_creator.py

├── shell_scripts

│   ├── indexing_reference_salmon.sh

│   └── quant_salmon.sh

└── README.md

## 📜 Script Descriptions
``` Root Directory ```
* ```Diff_Exp_Approach_II.R```: This is the central script of the entire project.

    * Inputs: Salmon quantification files (quant.sf) and a sample metadata table.
    
    * Actions:

        1. Imports transcript counts and aggregates them to the gene level using tximport.
        
        2. Creates DESeqDataSet objects for different experimental designs (stage, site, collection type).
        
        3. Runs the DESeq2 pipeline to identify differentially expressed genes.
        
        4. Saves the complete, shrunken results for all planned comparisons to CSV files.
        
        5. Performs variance stabilizing transformation (VST) and generates PCA plots to visualize sample relationships.
        
        6. Prepares ranked gene lists and runs Gene Set Enrichment Analysis (GSEA) for each comparison.

    * Outputs: Full and significant gene lists (CSV), PCA plots (PNG), and GSEA results (CSV and plots).

```/Reference_transcriptome```
This directory contains scripts for preparing custom annotation files required for the analysis.

  * tx2gene_creator.py: A Python script to parse a FASTA header and create a two-column CSV file mapping transcript accessions to gene IDs.
  
  * gaf2orgDB.R: An R script that reads a Gene Ontology Annotation File (GAF) and uses the AnnotationForge package to build a custom, installable OrgDb annotation package for Anopheles stephensi.

```/shell_scripts```
This directory contains shell scripts for running command-line tools, intended for a High-Performance Computing (HPC) environment using a PBS scheduler.

  * indexing_reference_salmon.sh: Creates a Salmon index from the reference transcriptome FASTA file. This only needs to be run once.
  
  * quant_salmon.sh: Iterates over a list of sample IDs, running salmon quant to quantify transcript abundance from paired-end FASTQ files.

```/Approach_II_outputs```
This directory contains scripts for downstream functional analysis that use the differential expression results.

  * ```GO_enrichment.R```: Performs over-representation analysis for Gene Ontology (GO) terms.

      * Inputs: Significant up- and down-regulated gene lists from Diff_Exp_Approach_II.R.
      
      * Actions: Uses clusterProfiler and the custom org.Astephensi.eg.db to find enriched GO terms for both unfiltered and directionally filtered gene sets.
      
      * Outputs: CSV files of enriched terms and corresponding bar plots/dot plots.

  * ```KEGGS_analysis.R```: Performs over-representation analysis for KEGG pathways.
  
      * Inputs: Significant up- and down-regulated gene lists.
      
      * Actions: Uses clusterProfiler's enrichKEGG function to find enriched pathways.
      
      * Outputs: CSV files of enriched pathways and corresponding plots.

```/Approach_II_outputs/Pathway_selection```
  * ```Pathway_selections.R```: A final analysis script to synthesize results.
  
      * Inputs: The enrichment results generated by GO_enrichment.R and KEGGS_analysis.R.
      
      * Actions: Loads the various pathway datasets into a nested list structure, filters out pathways associated with larval stage development, and then finds the intersection between pathways identified through different filtering strategies (gene-level vs. pathway-level filtering).
      
      * Outputs: CSV files containing the final, intersected lists of significant pathways for each comparison.

## 🛠️ Requirements
### Software
  * R (v4.1 or later)
  
  * Python (v3.6 or later)
  
  * Salmon (v1.9 or later)

  * Apptainer (formerly Singularity) (for containerized execution on HPC)

### R Packages
You can install the required R packages by running the following commands in an R session:

`R`
```
# CRAN Packages
install.packages(c("tidyverse", "readr", "dplyr", "tidyr", "stringr", "ggplot2", 
                   "pheatmap", "ggrepel", "patchwork", "RColorBrewer"))

# Bioconductor Packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("rtracklayer", "GenomicFeatures", "DESeq2", "tximport", 
                       "AnnotationDbi", "ashr", "clusterProfiler", "enrichplot",
                       "AnnotationForge", "txdbmaker", "biomaRt", "GO.db"))

# Note: The custom org.Astephensi.eg.db package is built and installed by the
# 'gaf2orgDB.R' script.

```
