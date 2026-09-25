setwd("/path/to/project/scripts/")

#Load the necessary libraries
library(Seurat)
library(qs2)

##====SUBSET TO SHARED GENE SPACE (the external controls and the internal samples have different gene sets-take intersect====
obj <- qs_read("../objects/obj_merge_filt_doublets_removed.qs2")
# shared symbols between internal  and external (GSE247531) references
base_dir  <- "/path/to/project/data/"
srsf2_sym <- read.table(gzfile(file.path(base_dir,"P6/filtered_feature_bc_matrix/features.tsv.gz")), sep="\t")$V2
ctrl_sym  <- read.table(gzfile(file.path(base_dir,"HD3/filtered_feature_bc_matrix/features.tsv.gz")), sep="\t")$V2
shared    <- intersect(unique(srsf2_sym), unique(ctrl_sym))
obj       <- subset(obj, features = shared)
cat("Genes after subsetting to shared space:", nrow(obj), "\n")

##====NORMALIZATION====
obj<- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
print("Data succesfully normalized")
qs_save(obj, "../objects/obj_merge_filt_doublets_removed_norm.qs2")

##====VARIABLE FEATURES====
obj <- qs_read("../objects/obj_merge_filt_doublets_removed_norm.qs2")
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)

##====SCALING====
#check RAM before running (better use more than 128 GB of RAM and delete unecessary objects from environment)
all.genes <- rownames(obj)
obj <- ScaleData(obj, features = all.genes)
rm("all.genes")
qs_save(obj, "../objects/obj_merge_filt_doublets_removed_norm_scaled.qs2")
print("Data succesfully scaled")
