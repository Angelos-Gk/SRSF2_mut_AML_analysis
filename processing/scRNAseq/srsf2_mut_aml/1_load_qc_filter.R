setwd("/path/to/project/scripts/")

#====LIBRARIES====
#Load the necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(SingleR)
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
library(qs2)

#====PREPROCESSING====
##====DATA LOADING====
#Define the base directory
#Expects one subdirectory per sample, named by sample ID (P1-P17, HD1-HD6),
#each containing the Cell Ranger output folder "filtered_feature_bc_matrix".
#HD3-HD6 are the healthy donor CD34+ samples HD1-HD4 from GEO accession GSE247531.
base_dir <- "/path/to/project/data/"
#List all subdirectories in the base directory
pool_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
#Initialize an empty list to store the Seurat objects
seurat_list <- list()
#Initialize the progress bar
pb <- progress_bar$new(
  format = "  Processing [:bar] :percent | Elapsed: :elapsedfull | ETA: :eta | Current: :current | Total: :total",
  total = length(pool_dirs), clear = FALSE
)
#Loop through each subdirectory, load the data, and create a Seurat object
for (i in seq_along(pool_dirs)) {
  dir <- pool_dirs[i]
  pool_name <- basename(dir)
  # Construct the path to the filtered_feature_bc_matrix directory
  data_dir <- file.path(dir, "filtered_feature_bc_matrix")
  #Load the data
  counts <- Read10X(data.dir = data_dir)
  if ("Gene Expression" %in% names(counts)) {
    counts <- counts$'Gene Expression'
  } else {
    counts <- counts
  }
  #Create a Seurat object
  seurat_obj <- CreateSeuratObject(counts = counts, project = pool_name, min.cells = 3)
  #Store the Seurat object in the list with the pool name as the key
  seurat_list[[pool_name]] <- seurat_obj
  #Update the progress bar
  pb$tick()
}
print("Data succesfully loaded")

##====MERGE SEURAT OBJECTS====
# Extract the names of all the Seurat objects
seurat_names <- names(seurat_list)

# Merge all Seurat objects in the list into one
obj <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],  # All other Seurat objects
  add.cell.ids = seurat_names,
  project = "srsf2"
)
#remove unecessary objects
rm(counts,pb,seurat_list,seurat_obj,base_dir,data_dir,dir,i,pool_dirs,pool_name,seurat_names)
#add info regarding mit content
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
#check info regarding the merged object
obj
#calculate how many cells from each sample
table(obj@meta.data$orig.ident)
#save the object
qs_save(obj, "../objects/obj_merge_raw.qs2")
cat("Saved with qs2\n")

##====QC====
obj <- qs_read("../objects/obj_merge_raw.qs2")
#Visualize QC metrics as a violin plot
Idents(obj) <- "orig.ident"
p1 <- SCpubr::do_ViolinPlot(sample = obj,
                            features = "percent.mt",
                            y_cut = 10,
                            legend.position = "none",
                            boxplot_width = 0.1)
p2 <- SCpubr::do_ViolinPlot(sample = obj,
                            features = "nFeature_RNA",
                            y_cut = 6000,
                            legend.position = "none",
                            boxplot_width = 0.1)
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
obj <- subset(obj,
              subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)

#check how many cells are filtered out from each sample
table(obj@meta.data$orig.ident)

print("Data succesfully filtered")

#====SAMPLE COLUMN====
# Sample folders are already named by sample ID, so the sample column equals orig.ident
obj@meta.data$sample <- obj@meta.data$orig.ident
head(obj@meta.data)
qs_save(obj, "../objects/obj_merge_filt.qs2")

#do doublet removal before continuing
