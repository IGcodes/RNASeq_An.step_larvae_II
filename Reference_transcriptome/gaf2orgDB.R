## Setting up working directory
setwd("C:/Users/aaisu/Box/Carter Lab/Projects/RNAseq_EthiopiaNUCI/Reference_transcriptome")

## --- 0) Setup ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("AnnotationForge","AnnotationDbi","GO.db"), ask = FALSE, update = FALSE)
install.packages(c("readr","dplyr","tidyr","stringr"))

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(AnnotationForge)
library(AnnotationDbi)
library(GO.db)

## --- 1) Inputs (EDIT THESE) ----
gaf_file <- "./UCI_ANSTEP_V1.0_gene_ontology.gaf"  # GAF 2.2
genus    <- "Anopheles"
species  <- "stephensi"
tax_id   <- "30069"
maint    <- "Isuru Gunarathna <isuru_gunarathna1@baylor.edu>" # Follow this exact format otherwise you're doomed!!!
author   <- "Isuru Gunarathna"
pkg_ver  <- "0.1.0"
outdir   <- "orgdb_from_gaf22_out"
drop_NOT <- TRUE   # set FALSE if you want to keep NOT annotations

dir.create(outdir, showWarnings = FALSE)

## --- 2) Read GAF 2.2 (skip comment lines starting with '!') ----
## Official GAF 2.2 columns (17)
gaf_cols <- c(
  "DB", "GeneID", "Symbol", "Qualifier", "GO_ID",
  "Reference", "Evidence_Code", "With_From", "Aspect",
  "Gene_Name", "Gene_Synonym", "Type", "Taxon",
  "Date", "Assigned_By", "Annot_Ext", "Gene_Product_Form_ID"
)

gaf <- readr::read_tsv(
  gaf_file,
  comment = "!",
  col_names = gaf_cols,
  col_types = cols(.default = "c"),
  na = c("", "NA"),
  progress = FALSE
) %>%
  filter(!is.na(GeneID), !is.na(GO_ID)) %>%
  mutate(
    GeneID       = as.character(GeneID),
    Symbol       = coalesce(Symbol, ""),
    Gene_Name    = coalesce(Gene_Name, ""),
    Type         = coalesce(Type, ""),
    Evidence_Code= coalesce(Evidence_Code, ""),
    Aspect       = coalesce(Aspect, ""),
    Qualifier    = coalesce(Qualifier, "")
  )

## Optionally drop NOT qualifiers
if (isTRUE(drop_NOT)) {
  gaf <- gaf %>% filter(!str_detect(Qualifier, "(^|\\|)NOT($|\\|)"))
}

## --- 3) Build gene_info (primary key MUST be named 'GID') ----
gene_info <- gaf %>%
  transmute(
    GID      = GeneID,
    SYMBOL   = if_else(Symbol == "" | is.na(Symbol), GeneID, Symbol),
    GENENAME = Gene_Name,
    GENTYPE  = Type
  ) %>%
  mutate(across(everything(), ~replace_na(.x, ""))) %>%
  distinct() %>%
  filter(GID != "")

## --- 4) Build GO table: GID, GO, EVIDENCE, ONTOLOGY ----
aspect_map <- c(P = "BP", C = "CC", F = "MF")
go_df <- gaf %>%
  transmute(
    GID      = GeneID,
    GO       = GO_ID,
    EVIDENCE = Evidence_Code,
    ONTOLOGY = unname(aspect_map[Aspect])
  ) %>%
  filter(GID != "", GO != "", EVIDENCE != "", !is.na(ONTOLOGY)) %>%
  distinct()

## --- 5) Sanity checks: column-name clashes (non-GID names must be unique across tables) ----
tables <- list(gene_info = gene_info, go = go_df)
non_gid_names <- lapply(tables, function(df) setdiff(colnames(df), "GID"))
dup_names <- Reduce(intersect, non_gid_names)
if (length(dup_names)) {
  stop("Column name(s) appear in multiple tables (other than GID): ",
       paste(dup_names, collapse = ", "),
       "\nRename them so each non-GID column is unique.")
}

## --- 6) Build the OrgDb package ----
args <- list(
  gene_info = gene_info,
  go        = go_df,
  version   = pkg_ver,
  maintainer= maint,
  author    = author,
  outputDir = outdir,
  tax_id    = tax_id,
  genus     = genus,
  species   = species
)

do.call(AnnotationForge::makeOrgPackage, args)
cat("Package written under: ", normalizePath(outdir), "\n")

## --- 7) Install & validate ----
pkg_dirs <- unique(dirname(list.files(outdir, pattern = "DESCRIPTION", recursive = TRUE, full.names = TRUE)))
if (length(pkg_dirs)) {
  install.packages(pkg_dirs, repos = NULL, type = "source")
}

## Load and check keytypes (GOALL should be present if GO.db was installed)
## Replace with your actual package name:
library(org.Astephensi.eg.db)
keytypes(org.Astephensi.eg.db)
columns(org.Astephensi.eg.db)
head(AnnotationDbi::select(org.Astephensi.eg.db,
                            keys = head(gene_info$GID, 3),
                            columns = c("SYMBOL","GO","EVIDENCE","ONTOLOGY"),
                            keytype = "GID"))
