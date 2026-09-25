setwd("/path/to/project/AML1026/scripts")

#====LIBRARIES====
#Load the necessary libraries
library(plyr)
library(Seurat)
library(ggplot2)
library(nanoranger.R)
library(dplyr)
library(patchwork)
library(SingleR)
library(celldex)
library(ComplexHeatmap)
library(png)
library(gridExtra)
library(grid)
library(tidyr)
library(SCpubr)
library(RColorBrewer)
library(SeuratWrappers)
library(progress)
library(reticulate)
library(sceasy)
library(devtools)
library(cellmarkeraccordion)

#Load the annotated AML1026 Seurat object
obj <- readRDS("../objects/obj_200_res1_final_anno.rds")

#====EXTRACT MUTATIONS====
# AML1026 nanorangerR extract - TWO DATASETS
#  (A) protocol_benchmark: only AML1026.1 with 3 protocols
#      - libraries: AML1026.1_rhPCR / AML1026.1_taqman / AML1026.1_taqman_biotin
#      - Seurat barcode in `bc` (AML1026.1_<CELLBC>-1)
#      - Unique row key in `bc_key` (bc__protocol)
#
#  (B) rhpcr_only: AML1026.1/2/3/4 rhPCR only
#      - libraries: AML1026.1_rhPCR / AML1026.2_rhPCR / AML1026.3_rhPCR / AML1026.4_rhPCR
#      - Seurat barcode in `bc` (AML1026.X_<CELLBC>-1)
#      - Unique row key in `bc_key` (bc__protocol) (protocol will be "rhPCR")
#
# Both datasets:
#  - skip header-only pileups (wc -l == 1)
#  - avoid duplicate rownames by using bc_key
#  - save per-gene TSVs + combined TSV + RDS list
#STARCODE location
Sys.setenv(PATH = paste(Sys.getenv("PATH"),
                        "/path/to/starcode",
                        sep=":"))

patient_id <- "AML1026"
base_out <- "/path/to/project/nanoranger/output/nanoranger"

#Define two library sets
libraries_benchmark <- c(
  "AML1026.1_rhPCR",
  "AML1026.1_taqman",
  "AML1026.1_taqman_biotin"
)

libraries_rhpcr <- c(
  "AML1026.1_rhPCR",
  "AML1026.2_rhPCR",
  "AML1026.3_rhPCR",
  "AML1026.4_rhPCR"
)

#Map each nanoranger folder -> Seurat prefix (must match Seurat rownames)
lib_to_seurat <- c(
  "AML1026.1_rhPCR"         = "AML1026.1",
  "AML1026.1_taqman"        = "AML1026.1",
  "AML1026.1_taqman_biotin" = "AML1026.1",
  "AML1026.2_rhPCR"         = "AML1026.2",
  "AML1026.3_rhPCR"         = "AML1026.3",
  "AML1026.4_rhPCR"         = "AML1026.4"
)

#Protocol is everything after the first underscore
get_protocol <- function(library_id) sub("^[^_]+_", "", library_id)

make_lib_tbl <- function(libraries) {
  data.frame(
    library_id    = libraries,
    seurat_prefix = unname(lib_to_seurat[libraries]),
    protocol      = vapply(libraries, get_protocol, character(1)),
    pileup_dir    = file.path(base_out, libraries, "pileups"),
    stringsAsFactors = FALSE
  )
}

lib_tbl_benchmark <- make_lib_tbl(libraries_benchmark)
lib_tbl_rhpcr      <- make_lib_tbl(libraries_rhpcr)

#Genes + ALT/REF
genes <- c("CCDC88B", "CSNK1D", "GNB1", "LPGAT1", "MARS1", "NRAS", "PNKP", "SF3B1", "SRSF2")

gene_info <- list(
  CCDC88B = c(ALT="C", REF="G"),
  CSNK1D  = c(ALT="C", REF="G"),
  GNB1    = c(ALT="C", REF="T"),
  LPGAT1  = c(ALT="T", REF="C"),
  MARS1   = c(ALT="T", REF="C"),
  NRAS    = c(ALT="G", REF="C"),
  PNKP    = c(ALT="C", REF="G"),
  SF3B1   = c(ALT="T", REF="C"),
  SRSF2   = c(ALT="T", REF="G")
)

#Output dirs 
out_base <- "/path/to/project/nanoranger/output/r_tables"
out_benchmark <- file.path(out_base, "protocol_benchmark")
out_rhpcr      <- file.path(out_base, "rhpcr_only")

dir.create(out_benchmark, recursive = TRUE, showWarnings = FALSE)
dir.create(out_rhpcr, recursive = TRUE, showWarnings = FALSE)

#Helpers 
pileup_path <- function(pileup_dir, gene) {
  file.path(pileup_dir, paste0(patient_id, "_pileup_", gene, ".csv.gz"))
}

#Header-only safe check: TRUE only if header + >=1 data row
has_data_rows_gz <- function(fp) {
  if (!file.exists(fp)) return(FALSE)
  nlines <- tryCatch({
    as.integer(system(paste0("zcat ", shQuote(fp), " | wc -l"), intern = TRUE))
  }, warning = function(w) NA_integer_,
  error   = function(e) NA_integer_)
  if (is.na(nlines)) return(FALSE)
  nlines > 1
}

safe_extract_mutation <- function(fp, ALT, REF, FILTER=4) {
  if (!has_data_rows_gz(fp)) return(NULL)
  tryCatch({
    df <- extract_mutation(fp, ALT, REF, FILTER = FILTER)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    as.data.frame(df)
  }, error = function(e) {
    message("[WARN] extract_mutation failed for: ", fp, " :: ", conditionMessage(e))
    NULL
  })
}

#Core runner: returns list(gene_list=..., mut_df=...) and writes per-gene TSVs 
run_extract_dataset <- function(lib_tbl, out_dir, dataset_name, FILTER=4) {
  mutation_data_list <- list()
  
  for (gene in genes) {
    gene_res <- list()
    ALT <- gene_info[[gene]][["ALT"]]
    REF <- gene_info[[gene]][["REF"]]
    
    for (i in seq_len(nrow(lib_tbl))) {
      lib <- lib_tbl$library_id[i]
      seurat_prefix <- lib_tbl$seurat_prefix[i]
      protocol <- lib_tbl$protocol[i]
      fp <- pileup_path(lib_tbl$pileup_dir[i], gene)
      
      df <- safe_extract_mutation(fp, ALT, REF, FILTER = FILTER)
      if (is.null(df)) next
      
      # Seurat-compatible barcode for matching to obj@meta.data rownames
      df$bc <- paste0(seurat_prefix, "_", df$bc, "-1")
      
      # Unique key to avoid duplicate rownames across protocols / libraries
      df$bc_key <- paste0(df$bc, "__", protocol)
      
      # Metadata
      df$dataset     <- dataset_name
      df$sample_name <- seurat_prefix
      df$protocol    <- protocol
      df$library_id  <- lib
      df$patient_id  <- patient_id
      df$gene        <- gene
      df$pileup_file <- fp
      
      gene_res[[lib]] <- df
    }
    
    if (length(gene_res) == 0) {
      message("[INFO] [", dataset_name, "] No usable pileup data for gene: ", gene)
      next
    }
    
    gene_df <- do.call(rbind, gene_res)
    gene_df <- as.data.frame(gene_df)
    rownames(gene_df) <- make.unique(gene_df$bc_key)
    
    write.table(gene_df,
                file = file.path(out_dir, paste0(patient_id, "_", gene, ".tsv")),
                quote = FALSE, sep = "\t", row.names = FALSE)
    
    mutation_data_list[[gene]] <- gene_df
  }
  
  if (length(mutation_data_list) == 0) {
    message("[WARN] [", dataset_name, "] No genes extracted at all.")
    mut_df <- data.frame()
  } else {
    mut_df <- do.call(rbind, mutation_data_list)
    mut_df <- as.data.frame(mut_df)
    rownames(mut_df) <- make.unique(mut_df$bc_key)
  }
  
  # Save list + combined
  saveRDS(mutation_data_list,
          file = file.path(out_dir, paste0(patient_id, "_mutation_data_list_", dataset_name, ".rds")))
  
  write.table(mut_df,
              file = file.path(out_dir, paste0(patient_id, "_ALL_GENES_", dataset_name, ".tsv")),
              quote = FALSE, sep = "\t", row.names = FALSE)
  
  message("[DONE] Dataset: ", dataset_name,
          " | genes: ", paste(names(mutation_data_list), collapse = ", "),
          " | rows: ", nrow(mut_df))
  
  list(gene_list = mutation_data_list, mut_df = mut_df)
}

#Run both datasets
res_benchmark <- run_extract_dataset(
  lib_tbl = lib_tbl_benchmark,
  out_dir = out_benchmark,
  dataset_name = "protocol_benchmark",
  FILTER = 4
)

res_rhpcr <- run_extract_dataset(
  lib_tbl = lib_tbl_rhpcr,
  out_dir = out_rhpcr,
  dataset_name = "rhpcr_only",
  FILTER = 4
)

#====MATCH WITH SEURAT====
save_dir <- "../objects"
stopifnot(dir.exists(save_dir))

# Seurat-valid barcodes
seurat_bcs <- rownames(obj@meta.data)

# --- helper: filter a nanoranger gene-list by Seurat barcodes ---
filter_gene_list_to_seurat <- function(gene_list, seurat_bcs, drop_empty = TRUE) {
  filtered <- lapply(gene_list, function(df) {
    df <- as.data.frame(df)
    if (!("bc" %in% colnames(df))) {
      stop("Expected column `bc` not found in one of the gene data.frames.")
    }
    df[df$bc %in% seurat_bcs, , drop = FALSE]
  })
  
  if (drop_empty) {
    filtered <- Filter(function(df) nrow(df) > 0, filtered)
  }
  filtered
}

# Load the per-dataset gene lists saved by run_extract_dataset() above
AML1026_benchmark <- readRDS(file.path(out_benchmark, "AML1026_mutation_data_list_protocol_benchmark.rds"))
AML1026_rhPCR     <- readRDS(file.path(out_rhpcr, "AML1026_mutation_data_list_rhpcr_only.rds"))

# Filter both datasets
AML1026_benchmark_filt <- filter_gene_list_to_seurat(AML1026_benchmark, seurat_bcs, drop_empty = TRUE)
AML1026_rhPCR_filt     <- filter_gene_list_to_seurat(AML1026_rhPCR,     seurat_bcs, drop_empty = TRUE)

# Save RDS
saveRDS(
  AML1026_benchmark_filt,
  file = file.path(save_dir, "AML1026_benchmark_mutation_data_list_inSeurat.rds")
)

saveRDS(
  AML1026_rhPCR_filt,
  file = file.path(save_dir, "AML1026_rhPCRonly_mutation_data_list_inSeurat.rds")
)
