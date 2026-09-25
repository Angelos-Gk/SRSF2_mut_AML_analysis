# Run AFTER: 01_load_qc_filter.R (only needs obj_merge_filt.qs2;
#            independent of 02_doubletfinder.R)
# Saves:     meta_scdblfinder.rds  (scDblFinder results per cell)
# Submit via: 03_scdblfinder_submit.sh
# ============================================================

# ---- Libraries ----
library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(dplyr)
library(qs2)

# ---- Paths ----
setwd("/path/to/project/scripts")
save_dir <- "../objects/"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load filtered object ----
obj <- qs_read(paste0(save_dir, "obj_merge_filt.qs2"))
cat("Loaded object with", ncol(obj), "cells\n")

# ---- Convert to SCE and run scDblFinder ----
# Coerce RNA assay to Assay (v3) for compatibility
obj[["RNA"]] <- as(obj[["RNA"]], Class = "Assay")
sce <- as.SingleCellExperiment(obj)

# Run scDblFinder per sample using the 'orig.ident' column
# dbr = NULL lets scDblFinder estimate the doublet rate automatically per sample
sce <- scDblFinder(sce, samples = "orig.ident", dbr = NULL)
cat("scDblFinder complete\n")
table(sce$scDblFinder.class)

# ---- Extract and save results ----
meta_scdblfinder <- as.data.frame(sce@colData) %>%
  dplyr::select(dplyr::starts_with("scDblFinder"))
rownames(meta_scdblfinder) <- sce@colData@rownames
head(meta_scdblfinder)

saveRDS(meta_scdblfinder, paste0(save_dir, "meta_scdblfinder.rds"))
cat("scDblFinder results saved to meta_scdblfinder.rds\n")
