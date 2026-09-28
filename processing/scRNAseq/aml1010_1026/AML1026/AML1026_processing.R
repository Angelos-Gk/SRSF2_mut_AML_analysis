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
library(symphony)
library(BoneMarrowMap)
library(qs2)

#====PREPROCESSING====
##====DATA LOADING====
#Read in each sample
AML1026.1 <- Read10X(data.dir = "../data/AML1026.1")
AML1026.2 <- Read10X(data.dir = "../data/AML1026.2")
AML1026.3 <- Read10X(data.dir = "../data/AML1026.3")
AML1026.4 <- Read10X(data.dir = "../data/AML1026.4")

#Create Seurat object for each sample (keep only gene expression assay)
AML1026.1 <- CreateSeuratObject(counts = AML1026.1$'Gene Expression', project = "AML1026.1", min.cells = 3)
AML1026.2 <- CreateSeuratObject(counts = AML1026.2$'Gene Expression', project = "AML1026.2", min.cells = 3)
AML1026.3 <- CreateSeuratObject(counts = AML1026.3$'Gene Expression', project = "AML1026.3", min.cells = 3)
AML1026.4 <- CreateSeuratObject(counts = AML1026.4$'Gene Expression', project = "AML1026.4", min.cells = 3)

#Merge the Seurat objects
obj <- merge(AML1026.1, y = c(AML1026.2, AML1026.3, AML1026.4), 
             add.cell.ids = c("AML1026.1","AML1026.2","AML1026.3","AML1026.4"), 
             project = "AML")

#remove unecessary objects
rm("AML1026.1", "AML1026.2", "AML1026.3", "AML1026.4")
#add info regarding mit content
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
#check info regarding the merged object
obj
#calculate how many cells from each sample
table(obj@meta.data$orig.ident)

# Rename using recode
obj@meta.data$orig.ident <- dplyr::recode(obj@meta.data$orig.ident,
                                   `AML1026.1` = "Screening",
                                   `AML1026.2` = "EOLN",
                                   `AML1026.3` = "C1",
                                   `AML1026.4` = "EOT")

# Verify the changes by viewing the updated table
table(obj@meta.data$orig.ident)

##====QC====
#Visualize QC metrics as a violin plot
Idents(obj) <- "orig.ident"
p1 <- SCpubr::do_ViolinPlot(sample = obj, 
                            features = "percent.mt", 
                            y_cut = 10,
                            legend.position = "none",
                            boxplot_width = 0.1)
p2 <- SCpubr::do_ViolinPlot(sample = obj, 
                            features = "nFeature_RNA", 
                            y_cut = 7000,  # This will apply the y_cut at 7000
                            legend.position = "none",
                            boxplot_width = 0.1)
# Add another horizontal line at y = 150, with custom dash length
p2 <- p2 + ggplot2::geom_hline(yintercept = 150, 
                               linetype = "longdash",  # Use a predefined or custom dash style
                               color = "black", size = 1.2)
p3 <- SCpubr::do_ViolinPlot(sample = obj, 
                            features = "nCount_RNA", 
                            legend.position = "none",
                            boxplot_width = 0.1)
plot <- p1 / p2 / p3
ggsave("../plots/QC/QC.png", plot, width = 14, height = 20, dpi = 300)
#remove unecessary objects
rm("p1","p2","p3","plot")

##====FILTERING====
#filter the object
obj <- subset(obj, subset = nFeature_RNA > 150 & 
                nFeature_RNA < 7000 & 
                percent.mt < 10)
#check how many cells are filtered out from each sample
table(obj@meta.data$orig.ident)

##====NORMALIZATION====
obj<- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

##====VARIABLE FEATURES====
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)

##====SCALING====
#check RAM before running (better use more than 128 GB of RAM and delete unecessary objects from environment)
all.genes <- rownames(obj)
obj <- ScaleData(obj, features = all.genes)
rm("all.genes")
saveRDS(obj, file = "../objects/obj_merge_scaled.rds")
##====INTEGRATION & CLUSTERING====
top2000 <- head(VariableFeatures(obj), 2000)
obj.2000 <- obj[top2000]
saveRDS(obj.2000, file = "../objects/obj.2000.rds")
use_condaenv("scvi-env", required = TRUE)
print("conda env good")
sc <- import("scanpy", convert = FALSE)
print("scanpy good")
scvi <- import("scvi", convert = FALSE)
print("scvi good")
jax <- import("jax", convert = FALSE)
print("jax good")
py_run_string("import sys; print(sys.version)")
Sys.setenv(JAX_PLATFORMS = "gpu")
obj.2000[["RNA"]] <- as(object = obj.2000[["RNA"]], Class = "Assay")
Data <- convertFormat(obj.2000, from = "seurat", to = "anndata", main_layer = "counts", drop_single_values = FALSE)
# Setup model
scvi$model$SCVI$setup_anndata(Data, batch_key="orig.ident")
print("setup good")
model <- scvi$model$SCVI(Data)
#  model
print("starting training good")
model$train(accelerator = 'gpu', devices = 'auto', max_epochs = as.integer(100))
print("training over good")
latent <- scvi$model$SCVI$get_latent_representation(model)
x <- as.data.frame(colnames(obj.2000))
latent_r <- py_to_r(latent)
latent_df <- as.data.frame(latent_r)
rownames(latent_df) <- x$`colnames(obj.2000)`
latent_matrix <- as.matrix(latent_df)
print("manipulations good")
obj[["scvi"]] <- CreateDimReducObject(embeddings = latent_matrix,
                                      key = "scvi_", assay = DefaultAssay(obj.2000))
print("embedding good")
obj <- FindNeighbors(obj, dims = 1:10, reduction = "scvi")
print("neighbors good")
obj <- FindClusters(obj, resolution = 2)
print("clusters good")
obj <- RunUMAP(obj, dims = 1:10, reduction = "scvi", n.components = 2)
print("Data succesfully integrated")

saveRDS(obj, file = "../objects/obj_merge_scvi.rds")
print("save object: Good")

##====CELL TYPE ANNOTATION (BoneMarrowMap)====
# -------------------------
# Reference + uwot model
# -------------------------
# Folder with the BoneMarrowMap reference files (BoneMarrowMap_SymphonyReference.rds and
# BoneMarrowMap_uwot_model.uwot), downloaded as described in the BoneMarrowMap README
projection_path <- "/path/to/BoneMarrowMap_reference"

ref <- readRDS(file.path(projection_path, "BoneMarrowMap_SymphonyReference.rds"))
ref$save_uwot_path <- file.path(projection_path, "BoneMarrowMap_uwot_model.uwot")

# -------------------------
# Helper: run BMMap on one Seurat object and ensure ALL cells get labels
# -------------------------
annotate_with_bonemarrowmap_all_cells <- function(
  obj,
  batchvar = NULL,                       # e.g. "sample" if present; NULL if none
  mad_threshold = 10,                   # keep QC metric; can increase to be more permissive
  final_label = "BoneMarrowMap_CellType" # column name you want
) {

  # raw counts
  exp_mat <- Seurat::GetAssayData(obj, assay = "RNA", layer = "counts")

  # map query onto reference
  obj_mapped <- BoneMarrowMap::map_Query(
    exp_query      = exp_mat,
    metadata_query = obj@meta.data,
    ref_obj        = ref,
    vars           = batchvar
  )

  # calculate mapping error QC (does NOT remove cells; just annotates pass/fail)
  obj_mapped <- obj_mapped %>%
    BoneMarrowMap::calculate_MappingError(reference = ref, MAD_threshold = mad_threshold)

  # predict cell types:
  # - initial_label: assigned to ALL cells (pre-QC)
  # - final_label: will be NA for QC-fail cells
  obj_mapped <- BoneMarrowMap::predict_CellTypes(
    query_obj    = obj_mapped,
    ref_obj      = ref,
    initial_label = paste0(final_label, "_initial"),
    final_label   = final_label
  )

  # FORCE a label for every cell:
  # if final label is NA (QC fail), replace with initial label
  obj_mapped@meta.data[[final_label]] <- dplyr::coalesce(
    obj_mapped@meta.data[[final_label]],
    obj_mapped@meta.data[[paste0(final_label, "_initial")]]
  )

  # optional: also create a single "Broad" column filled everywhere
  # (BoneMarrowMap usually creates *_Broad and *_Broad_initial too)
  broad_final <- paste0(final_label, "_Broad")
  broad_init  <- paste0(final_label, "_Broad_initial")
  if (broad_final %in% colnames(obj_mapped@meta.data) && broad_init %in% colnames(obj_mapped@meta.data)) {
    obj_mapped@meta.data[[broad_final]] <- dplyr::coalesce(
      obj_mapped@meta.data[[broad_final]],
      obj_mapped@meta.data[[broad_init]]
    )
  }

  obj_mapped
}

# -------------------------
# Run for AML1026
# -------------------------
# If you have a batch variable in meta.data, set it here (e.g. "sample").
# If not, leave as NULL.
batchvar <- NULL

# join the per-sample count layers created by merge() so the full counts matrix can be mapped
obj <- JoinLayers(obj)

obj_1026_bmmap <- annotate_with_bonemarrowmap_all_cells(
  obj = obj,
  batchvar = batchvar,
  mad_threshold = 10,
  final_label = "BoneMarrowMap_CellType"
)

# -------------------------
# Quick sanity checks: are there any NAs left?
# -------------------------
table(is.na(obj_1026_bmmap$BoneMarrowMap_CellType))
table(is.na(obj_1026_bmmap$BoneMarrowMap_CellType_Broad))

# QC distribution (so you still know which cells were low-confidence)
table(obj_1026_bmmap$mapping_error_QC, useNA = "ifany")
table(obj_1026_bmmap$BoneMarrowMap_CellType, useNA = "ifany") |> sort(decreasing = TRUE)
table(obj_1026_bmmap$BoneMarrowMap_CellType_Broad, useNA = "ifany") |> sort(decreasing = TRUE)

# 1) Define the mapping function for the broad final annotation
make_bmmap_broad_final <- function(broad_vec) {
  dplyr::case_when(
    broad_vec %in% c("HSC MPP", "Cycling Progenitor") ~ "HSCs",
    broad_vec %in% c("Early GMP", "Late GMP") ~ "GMPs",
    broad_vec %in% c("Early Lymphoid", "LMPP") ~ "Lymphoid progenitors",
    broad_vec %in% c("Early Erythroid", "Late Erythroid") ~ "Erythroid-like cells",
    broad_vec %in% c("Megakaryocyte Precursor", "MEP", "EoBasoMast Precursor") ~ "Erythroid progenitors",
    broad_vec %in% c("B", "Pre-B", "Pro-B") ~ "B-cells",
    broad_vec %in% c("Naive T", "CD4 Memory T", "CD8 Memory T") ~ "T-cells",
    broad_vec %in% c("NK") ~ "NK cells",
    broad_vec %in% c("Pro-Monocyte", "Monocyte") ~ "Monocytes",
    broad_vec %in% c("cDC", "pDC") ~ "Dendritic cells",
    broad_vec %in% c("Stromal") ~ "Stromal cells",
    broad_vec %in% c("Plasma Cell") ~ "Plasma cells",
    TRUE ~ as.character(broad_vec)
  )
}

# 2) Apply to AML1026
obj_1026_bmmap$BoneMarrowMap_CellType_Broad_final <-
  make_bmmap_broad_final(obj_1026_bmmap$BoneMarrowMap_CellType_Broad)

# 3) Fix the level order (important for downstream comparisons/plots)
FINAL_LEVELS <- c(
  "HSCs",
  "GMPs",
  "Lymphoid progenitors",
  "Erythroid-like cells",
  "Erythroid progenitors",
  "Monocytes",
  "T-cells",
  "B-cells",
  "NK cells",
  "Dendritic cells",
  "Plasma cells",
  "Stromal cells"
)

obj_1026_bmmap$BoneMarrowMap_CellType_Broad_final <-
  factor(obj_1026_bmmap$BoneMarrowMap_CellType_Broad_final, levels = FINAL_LEVELS)

# 4) Sanity checks
cat("\nAML1026 Broad_final counts:\n")
print(table(obj_1026_bmmap$BoneMarrowMap_CellType_Broad_final, useNA = "ifany"))

cat("\nAny categories missing in AML1026? (0 counts)\n")
print(FINAL_LEVELS[table(obj_1026_bmmap$BoneMarrowMap_CellType_Broad_final) == 0])

qs2::qs_save(obj_1026_bmmap, "../objects/obj_1026_bonemarrowmap_allcells.qs2")

#transfer annotation to the main object
# 1) Sanity check: do the cell names match?
stopifnot(all(colnames(obj) %in% colnames(obj_1026_bmmap)))

# 2) Transfer the column (name-aligned)
obj$BoneMarrowMap_CellType_Broad_final <-
  obj_1026_bmmap$BoneMarrowMap_CellType_Broad_final[colnames(obj)]

# 3) Quick checks
table(is.na(obj$BoneMarrowMap_CellType_Broad_final))
table(obj$BoneMarrowMap_CellType_Broad_final, useNA = "ifany")

qs2::qs_save(obj, file = "../objects/obj_1026_with_BoneMarrowMap_Broad_final.qs2")
