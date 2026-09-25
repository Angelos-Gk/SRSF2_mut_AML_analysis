# Run AFTER: 01_load_qc_filter.R (obj_merge_filt.qs2)
# Saves:     sglt_dblt_metadata_DoubletFinder.rds  (DoubletFinder results per cell)
# ============================================================

# ---- Libraries ----
library(Seurat)
library(DoubletFinder)
library(dplyr)
library(qs2)

options(max.print = .Machine$integer.max, scipen = 999,
        stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE)

# ---- Paths ----
setwd("/path/to/project/scripts")
save_dir    <- "../objects/"
plot_dir    <- "../plots/Doublet_detection/"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load filtered object ----
obj <- qs_read(paste0(save_dir, "obj_merge_filt.qs2"))
obj$SampleID <- obj$orig.ident
Idents(obj) <- "SampleID"
cat("Loaded object with", ncol(obj), "cells across", length(unique(obj$orig.ident)), "samples\n")

# ---- DoubletFinder custom function ----
run_doubletfinder_custom <- function(seu_sample_subset, multiplet_rate = NULL) {

  sample_name <- unique(seu_sample_subset$orig.ident)
  print(paste0("Sample ", sample_name, "..........."))

  # Estimate multiplet rate from 10x table if not provided (base R to avoid dplyr conflicts)
  if (is.null(multiplet_rate)) {
    multiplet_rates_10x <- data.frame(
      Multiplet_rate  = c(0.004, 0.008, 0.016, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
      Recovered_cells = c(500,   1000,  2000,  3000,  4000,  5000,  6000,  7000,  8000,  9000,  10000)
    )
    n_cells     <- nrow(seu_sample_subset@meta.data)
    valid_rows  <- multiplet_rates_10x[multiplet_rates_10x$Recovered_cells <= n_cells, ]
    if (nrow(valid_rows) == 0) {
      multiplet_rate <- 0.004
    } else {
      multiplet_rate <- valid_rows$Multiplet_rate[which.max(valid_rows$Recovered_cells)]
    }
    print(paste("Setting multiplet rate to", multiplet_rate, "for", n_cells, "cells"))
  }

  # Pre-process for DoubletFinder
  sample <- NormalizeData(seu_sample_subset, verbose = FALSE)
  sample <- FindVariableFeatures(sample, verbose = FALSE)
  sample <- ScaleData(sample, verbose = FALSE)
  sample <- RunPCA(sample, nfeatures.print = 0, verbose = FALSE)

  # Automatically select significant PCs
  stdv        <- sample[["pca"]]@stdev
  percent_stdv <- (stdv / sum(stdv)) * 100
  cumulative  <- cumsum(percent_stdv)
  co1 <- which(cumulative > 90 & percent_stdv < 5)[1]
  co2 <- sort(which((percent_stdv[1:(length(percent_stdv) - 1)] -
                       percent_stdv[2:length(percent_stdv)]) > 0.1),
              decreasing = TRUE)[1] + 1
  min_pc <- min(co1, co2)
  if (is.na(min_pc) || min_pc < 5) min_pc <- 10  # fallback
  print(paste("Using", min_pc, "PCs"))

  # Finish preprocessing
  sample <- RunUMAP(sample, dims = 1:min_pc, verbose = FALSE)
  sample <- FindNeighbors(sample, dims = 1:min_pc, verbose = FALSE)
  sample <- FindClusters(sample, resolution = 0.1, verbose = FALSE)

  # pK identification
  sweep_list  <- paramSweep(sample, PCs = 1:min_pc, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn       <- find.pK(sweep_stats)
  optimal.pk  <- as.numeric(as.character(
    bcmvn$pK[which.max(bcmvn$BCmetric)]
  ))

  # Homotypic doublet proportion estimate
  annotations    <- sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp.poi       <- round(multiplet_rate * nrow(sample@meta.data))
  nExp.poi.adj   <- round(nExp.poi * (1 - homotypic.prop))

  # Run DoubletFinder
  sample <- doubletFinder(seu  = sample,
                          PCs  = 1:min_pc,
                          pK   = optimal.pk,
                          nExp = nExp.poi.adj)

  # Standardise column name
  colnames(sample@meta.data)[grepl("DF.classifications", colnames(sample@meta.data))] <- "doublet_finder"

  # Return just the result dataframe
  result <- sample@meta.data["doublet_finder"]
  result <- tibble::rownames_to_column(result, "row_names")
  return(result)
}

# ---- Split and run per sample ----
samp_split <- SplitObject(obj, split.by = "SampleID")

samp_split_results <- lapply(samp_split, function(s) {
  tryCatch(
    run_doubletfinder_custom(s),
    error = function(e) {
      message("ERROR on sample ", unique(s$orig.ident), ": ", e$message)
      # Return NAs so the pipeline doesn't break
      result <- data.frame(row_names = colnames(s), doublet_finder = NA_character_)
      return(result)
    }
  )
})

# ---- Merge results and save ----
sglt_dblt_metadata <- data.frame(dplyr::bind_rows(samp_split_results))
rownames(sglt_dblt_metadata) <- sglt_dblt_metadata$row_names
sglt_dblt_metadata$row_names <- NULL
head(sglt_dblt_metadata)

saveRDS(sglt_dblt_metadata, paste0(save_dir, "sglt_dblt_metadata_DoubletFinder.rds"))
cat("DoubletFinder done. Results saved.\n")
sglt_dblt_metadata <- readRDS(paste0(save_dir, "sglt_dblt_metadata_DoubletFinder.rds"))
table(sglt_dblt_metadata$doublet_finder)

# ---- Summary per sample ----
obj <- AddMetaData(obj, sglt_dblt_metadata, col.name = "doublet_finder")

doublets_summary_df <- obj@meta.data %>%
  group_by(orig.ident, doublet_finder) %>%
  summarise(total_count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(countT   = sum(total_count),
         percent  = paste0(round(100 * total_count / countT, 2), "%")) %>%
  dplyr::select(-countT)

write.table(doublets_summary_df,
            file      = paste0(plot_dir, "DoubletFinder_summary.txt"),
            quote     = FALSE, row.names = FALSE, sep = "\t")
print(doublets_summary_df, n = Inf)

# ---- Visualise ----
library(ggplot2)
library(patchwork)

# Extract data manually
df <- FetchData(obj, vars = c("nFeature_RNA", "percent.mt",
                              "sample", "doublet_finder"),
                layer = "counts")

vln_plot <- function(yvar) {
  ggplot(df, aes(x = sample, y = .data[[yvar]], fill = doublet_finder)) +
    geom_violin(scale = "width", position = position_dodge(0.9)) +
    theme_classic() +
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = yvar)
}

p_vln <- vln_plot("nFeature_RNA") / vln_plot("percent.mt")
p_vln
