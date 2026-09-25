setwd("/path/to/project/AML1010/scripts")

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

#Load the annotated AML1010 Seurat object
obj <- readRDS("../objects/obj_200_res1_final_anno.rds")

#====EXTRACT MUTATION====
# STARCODE location
Sys.setenv(PATH = paste(Sys.getenv("PATH"),
                        "/path/to/starcode",
                        sep=":"))

#Patient-specific paths
patient_id <- "AML1010"
base_out <- "/path/to/project/nanoranger/output/nanoranger"

#Define two library sets
# (A) protocol benchmark: AML1010.1 only, with 3 protocols
libraries_benchmark <- c(
  "AML1010.1_rhPCR",
  "AML1010.1_taqman",
  "AML1010.1_taqman_biotin"
)

# (B) rhPCR only across AML1010.1/3/4/5
libraries_rhpcr <- c(
  "AML1010.1_rhPCR",
  "AML1010.3_rhPCR",
  "AML1010.4_rhPCR",
  "AML1010.5_rhPCR"
)

#Map each nanoranger folder -> Seurat prefix (must match Seurat rownames)
lib_to_seurat <- c(
  "AML1010.1_rhPCR"         = "AML1010.1",
  "AML1010.1_taqman"        = "AML1010.1",
  "AML1010.1_taqman_biotin" = "AML1010.1",
  "AML1010.3_rhPCR"         = "AML1010.3",
  "AML1010.4_rhPCR"         = "AML1010.4",
  "AML1010.5_rhPCR"         = "AML1010.5"
)

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
lib_tbl_rhpcr     <- make_lib_tbl(libraries_rhpcr)

#Genes + ALT/REF (AML1010)
snv_genes <- c("NRAS", "RNPEPL1", "C4orf33", "CYFIP2", "E2F3", "MCM7",
               "BMPR1A", "ACIN1", "RALY", "THAP7", "SPNS1", "ZNF407")

snv_info <- list(
  NRAS    = c(ALT="T", REF="C"),
  RNPEPL1 = c(ALT="C", REF="G"),
  C4orf33 = c(ALT="A", REF="G"),
  CYFIP2  = c(ALT="T", REF="C"),
  E2F3    = c(ALT="G", REF="A"),
  MCM7    = c(ALT="C", REF="A"),
  BMPR1A  = c(ALT="C", REF="A"),
  ACIN1   = c(ALT="A", REF="G"),
  RALY    = c(ALT="A", REF="G"),
  THAP7   = c(ALT="T", REF="C"),
  SPNS1   = c(ALT="C", REF="A"),
  ZNF407  = c(ALT="G", REF="T")
)

#Indels from your mutations.1010.csv
indel_genes <- c("TCF12", "TCF12_2")
indel_info <- list(
  TCF12   = list(REF="G", CONSENSUS=+1),
  TCF12_2 = list(REF="C", CONSENSUS=+1)
)

#Output dirs
out_base <- "/path/to/project/nanoranger/output/r_tables"
out_benchmark <- file.path(out_base, "protocol_benchmark")
out_rhpcr     <- file.path(out_base, "rhpcr_only")

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

safe_extract_indel <- function(fp, REF, CONSENSUS=+1) {
  if (!has_data_rows_gz(fp)) return(NULL)
  tryCatch({
    df <- extract_indel(fp, REF = REF, CONSENSUS = CONSENSUS)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    as.data.frame(df)
  }, error = function(e) {
    message("[WARN] extract_indel failed for: ", fp, " :: ", conditionMessage(e))
    NULL
  })
}

#Core runner: extracts SNVs + INDELs; returns list(gene_list=..., mut_df=...)
run_extract_dataset <- function(lib_tbl, out_dir, dataset_name, FILTER=4) {
  mutation_data_list <- list()
  
  #SNVs
  for (gene in snv_genes) {
    gene_res <- list()
    ALT <- snv_info[[gene]][["ALT"]]
    REF <- snv_info[[gene]][["REF"]]
    
    for (i in seq_len(nrow(lib_tbl))) {
      lib <- lib_tbl$library_id[i]
      seurat_prefix <- lib_tbl$seurat_prefix[i]
      protocol <- lib_tbl$protocol[i]
      fp <- pileup_path(lib_tbl$pileup_dir[i], gene)
      
      df <- safe_extract_mutation(fp, ALT, REF, FILTER = FILTER)
      if (is.null(df)) next
      
      df$bc <- paste0(seurat_prefix, "_", df$bc, "-1")
      df$bc_key <- paste0(df$bc, "__", protocol)
      
      df$dataset     <- dataset_name
      df$sample_name <- seurat_prefix
      df$protocol    <- protocol
      df$library_id  <- lib
      df$patient_id  <- patient_id
      df$gene        <- gene
      df$variant_type <- "SNV"
      df$pileup_file <- fp
      
      gene_res[[lib]] <- df
    }
    
    if (length(gene_res) == 0) {
      message("[INFO] [", dataset_name, "] No usable SNV pileup data for gene: ", gene)
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
  
  #INDELs
  for (gene in indel_genes) {
    gene_res <- list()
    REF <- indel_info[[gene]]$REF
    CONS <- indel_info[[gene]]$CONSENSUS
    
    for (i in seq_len(nrow(lib_tbl))) {
      lib <- lib_tbl$library_id[i]
      seurat_prefix <- lib_tbl$seurat_prefix[i]
      protocol <- lib_tbl$protocol[i]
      fp <- pileup_path(lib_tbl$pileup_dir[i], gene)
      
      df <- safe_extract_indel(fp, REF = REF, CONSENSUS = CONS)
      if (is.null(df)) next
      
      df$bc <- paste0(seurat_prefix, "_", df$bc, "-1")
      df$bc_key <- paste0(df$bc, "__", protocol)
      
      df$dataset     <- dataset_name
      df$sample_name <- seurat_prefix
      df$protocol    <- protocol
      df$library_id  <- lib
      df$patient_id  <- patient_id
      df$gene        <- gene
      df$variant_type <- "INDEL"
      df$pileup_file <- fp
      
      gene_res[[lib]] <- df
    }
    
    if (length(gene_res) == 0) {
      message("[INFO] [", dataset_name, "] No usable INDEL pileup data for gene: ", gene)
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
  
  #Combined
  if (length(mutation_data_list) == 0) {
    message("[WARN] [", dataset_name, "] No genes extracted at all.")
    mut_df <- data.frame()
  } else {
    mut_df <- do.call(rbind, mutation_data_list)
    mut_df <- as.data.frame(mut_df)
    rownames(mut_df) <- make.unique(mut_df$bc_key)
  }
  
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

# MATCH WITH SEURAT (AML1010) + SAVE FILTERED RDS
save_dir <- "../objects"
stopifnot(dir.exists(save_dir))

seurat_bcs <- rownames(obj@meta.data)

filter_gene_list_to_seurat <- function(gene_list, seurat_bcs, drop_empty = TRUE) {
  filtered <- lapply(gene_list, function(df) {
    df <- as.data.frame(df)
    if (!("bc" %in% colnames(df))) stop("Expected column `bc` not found.")
    df[df$bc %in% seurat_bcs, , drop = FALSE]
  })
  if (drop_empty) filtered <- Filter(function(df) nrow(df) > 0, filtered)
  filtered
}

AML1010_benchmark <- readRDS(file.path(out_benchmark, "AML1010_mutation_data_list_protocol_benchmark.rds"))
AML1010_rhPCR     <- readRDS(file.path(out_rhpcr, "AML1010_mutation_data_list_rhpcr_only.rds"))

AML1010_benchmark_filt <- filter_gene_list_to_seurat(AML1010_benchmark, seurat_bcs, drop_empty = TRUE)
AML1010_rhPCR_filt     <- filter_gene_list_to_seurat(AML1010_rhPCR,     seurat_bcs, drop_empty = TRUE)

saveRDS(AML1010_benchmark_filt,
        file = file.path(save_dir, "AML1010_benchmark_mutation_data_list_inSeurat.rds"))

saveRDS(AML1010_rhPCR_filt,
        file = file.path(save_dir, "AML1010_rhPCRonly_mutation_data_list_inSeurat.rds"))
