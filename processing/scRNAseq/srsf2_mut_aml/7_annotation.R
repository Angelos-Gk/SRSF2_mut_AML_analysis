setwd("/path/to/project/scripts/")

library(Seurat)
library(symphony)
library(dplyr)
library(BoneMarrowMap)
library(qs2)
library(ggplot2)

obj <- qs_read("../objects/obj_merge_filt_doublets_removed_norm_scaled_scvi.qs2")
print("Load:Good")
# Folder with the BoneMarrowMap reference files (BoneMarrowMap_SymphonyReference.rds and
# BoneMarrowMap_uwot_model.uwot), downloaded as described in the BoneMarrowMap README
projection_path <- "/path/to/BoneMarrowMap_reference/"
# Load reference and UMAP model
ref <- readRDS(file.path(projection_path, "BoneMarrowMap_SymphonyReference.rds"))
ref$save_uwot_path <- file.path(projection_path, "BoneMarrowMap_uwot_model.uwot")
obj <- JoinLayers(obj)
gc()
print("Join Layers:Good")
exp_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
print("Expression matrix:Good")
dim(exp_mat)   
head(rownames(exp_mat))
head(colnames(exp_mat))

batchvar <- "sample"   # or set to NULL if you don't want batch correction

obj <- map_Query(
  exp_query      = exp_mat,         # your counts matrix
  metadata_query = obj@meta.data,   # your cell metadata
  ref_obj        = ref,
  vars           = batchvar
)
print("Map Query:Good")
obj <- obj %>%
  calculate_MappingError(
    reference     = ref,
    MAD_threshold = 3
  )
print("Calculate Mapping Error:Good")

# Optional visualization
plot_MappingErrorQC(obj)
plots <- plot_MappingErrorQC(obj)

# Example: save individual plots as PNG
ggsave("../plots/MapQC_Histogram.png",  plot = plots$MapQC_Histogram,  width = 6, height = 4, dpi = 300)
ggsave("../plots/MapQC_RNAcounts.png",  plot = plots$MapQC_RNAcounts,  width = 6, height = 4, dpi = 300)
ggsave("../plots/MapQC_UMAP.png",       plot = plots$MapQC_UMAP,       width = 6, height = 4, dpi = 300)
ggsave("../plots/MapQC_Hemoglobin.png", plot = plots$MapQC_Hemoglobin, width = 6, height = 4, dpi = 300)
print("PLot Mapping Error:Good")

obj <- predict_CellTypes(
  query_obj   = obj,
  ref_obj     = ref,
  initial_label = "initial_CellType",
  final_label = "BoneMarrowMap_CellType"   
)
print("Annotate:Good")

#table(obj$BoneMarrowMap_CellType_Broad)
#table(is.na(obj$BoneMarrowMap_CellType_Broad))
#table(obj_pass$initial_CellType)
#table(obj_pass$BoneMarrowMap_CellType)
#table(obj_pass$BoneMarrowMap_CellType_Broad)

qs_save(obj, "../objects/ext_obj_merge_filt_doublets_removed_norm_scaled_scvi_bonemarmap.qs2")
print("Save:Good")

#====BoneMarrowMap-Custom Anno====
obj$BoneMarrowMap_CellType_Broad_final <- case_when(
  obj$BoneMarrowMap_CellType_Broad %in% c("HSC MPP", "Cycling Progenitor") ~ "HSCs",
  obj$BoneMarrowMap_CellType_Broad %in% c("Early GMP", "Late GMP","EoBasoMast Precursor" ) ~ "GMPs",
  obj$BoneMarrowMap_CellType_Broad %in% c("Early Lymphoid", "LMPP") ~ "Lymphoid progenitors",
  obj$BoneMarrowMap_CellType_Broad %in% c("Early Erythroid", "Late Erythroid") ~ "Erythroid-like cells",
  obj$BoneMarrowMap_CellType_Broad %in% c("Megakaryocyte Precursor", "MEP" ) ~ "MEPs",
  obj$BoneMarrowMap_CellType_Broad %in% c("B", "Pre-B", "Pro-B") ~ "B-cells",
  obj$BoneMarrowMap_CellType_Broad %in% c("Naive T", "CD4 Memory T", "CD8 Memory T" ) ~ "T-cells",
  obj$BoneMarrowMap_CellType_Broad %in% c("NK") ~ "NK cells",
  obj$BoneMarrowMap_CellType_Broad %in% c("Pro-Monocyte", "Monocyte") ~ "Monocytes",
  obj$BoneMarrowMap_CellType_Broad %in% c("cDC", "pDC") ~ "Dendritic cells",
  obj$BoneMarrowMap_CellType_Broad %in% c("Stromal") ~ "Stromal cells",
  obj$BoneMarrowMap_CellType_Broad %in% c("Plasma Cell") ~ "Plasma cells", 
  TRUE ~ as.character(obj$BoneMarrowMap_CellType_Broad)
)
print("Manual:Good")
table(obj$BoneMarrowMap_CellType_Broad_final)
qs_save(obj, "../objects/ext_obj_merge_filt_doublets_removed_norm_scaled_scvi_bonemarmap_finalanno.qs2")
print("Save Manual:Good")

#Label transfer
# #transfer the labels to the original object:
obj1 <- qs_read("../objects/obj_merge_filt_doublets_removed_norm_scaled_scvi.qs2")
obj2 <- qs_read("../objects/ext_obj_merge_filt_doublets_removed_norm_scaled_scvi_bonemarmap_finalanno.qs2")

## ---- 1. Keep only Pass cells in obj2 ----
cat("obj2 before:", ncol(obj2), "cells\n")
obj2 <- subset(obj2, subset = mapping_error_QC == "Pass")
cat("obj2 after Pass filter:", ncol(obj2), "cells\n")  

## ---- 2. Subset obj1 to the SAME Pass cells (align by barcode) ----
pass_cells <- colnames(obj2)
# safety: every Pass cell must exist in obj1
stopifnot(all(pass_cells %in% colnames(obj1)))
obj1 <- subset(obj1, cells = pass_cells)
cat("obj1 after subset:", ncol(obj1), "cells\n")         

## ---- 3. Transfer the 4 label columns, aligned by barcode ----
cols_to_add <- c("initial_CellType",
                 "BoneMarrowMap_CellType_Broad",
                 "BoneMarrowMap_CellType_No_Abbreviations",
                 "BoneMarrowMap_CellType_Broad_final")

# confirm all 4 exist in obj2
stopifnot(all(cols_to_add %in% colnames(obj2@meta.data)))

# build a data.frame indexed by cell barcode, reordered to obj1's cell order
meta_to_add <- obj2@meta.data[colnames(obj1), cols_to_add, drop = FALSE]
# alignment guarantee
stopifnot(identical(rownames(meta_to_add), colnames(obj1)))

obj1 <- AddMetaData(obj1, metadata = meta_to_add)

## ---- 4. Verify transfer ----
cat("\nColumns now in obj1:\n"); print(colnames(obj1@meta.data))
cat("\nNA check (should all be 0):\n")
print(colSums(is.na(obj1@meta.data[, cols_to_add])))
cat("\nCross-check obj1 vs obj2 broad_final tables (must match):\n")
print(identical(table(obj1$BoneMarrowMap_CellType_Broad_final),
                table(obj2$BoneMarrowMap_CellType_Broad_final)))

##---- 5. transfer umap_projected reduction ----
red_name <- "umap_projected"
stopifnot(red_name %in% names(obj2@reductions))
## 1. pull embeddings, reorder rows to obj1's cell order
emb <- Embeddings(obj2, reduction = red_name)
emb <- emb[colnames(obj1), , drop = FALSE]
## 2. alignment guarantee (same check you used for metadata)
stopifnot(identical(rownames(emb), colnames(obj1)))
## 3. preserve obj2's original key so plots label axes correctly
orig_key <- Key(obj2[[red_name]])          # e.g. "umapprojected_"
## 4. attach to obj1
obj1[[red_name]] <- CreateDimReducObject(
  embeddings = emb,
  key        = orig_key,
  assay      = DefaultAssay(obj1)
)
## ---- verify ----
cat("\nReductions now in obj1:\n"); print(names(obj1@reductions))
cat("\nProjected-UMAP dims:", paste(dim(Embeddings(obj1, red_name)), collapse = " x "), "\n")
stopifnot(identical(rownames(Embeddings(obj1, red_name)), colnames(obj1)))
cat("Row alignment OK\n")

## ---- save ----
qs_save(obj1, "../objects/ext_obj_merge_filt_doublets_removed_norm_scaled_scvi_bonemarmap_passonly_finalanno.qs2")
