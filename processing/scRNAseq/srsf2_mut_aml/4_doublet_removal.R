# Run AFTER: 02_doubletfinder.R and 03_scdblfinder.R have both completed
# Requires:  sglt_dblt_metadata_DoubletFinder.rds
#            meta_scdblfinder.rds
#            obj_merge_filt.qs2
# Saves:     obj_merge_filt_doublets_removed.qs2
# ============================================================

# ---- Libraries ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(qs2)

options(max.print = .Machine$integer.max, scipen = 999,
        stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE)

# ---- Paths ----
setwd("/path/to/project/scripts")
save_dir <- "../objects/"
plot_dir <- "../plots/Doublet_detection/"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load object and doublet results ----
obj <- qs_read(paste0(save_dir, "obj_merge_filt.qs2"))

sglt_dblt_metadata <- readRDS(paste0(save_dir, "sglt_dblt_metadata_DoubletFinder.rds"))
meta_scdblfinder   <- readRDS(paste0(save_dir, "meta_scdblfinder.rds"))

# ---- Add both to object ----
obj <- AddMetaData(obj, sglt_dblt_metadata, col.name = "doublet_finder")
obj <- AddMetaData(obj, meta_scdblfinder %>% dplyr::select("scDblFinder.class"))

head(obj@meta.data)
cat("doublet_finder table:\n");   print(table(obj$doublet_finder))
cat("scDblFinder.class table:\n"); print(table(obj$scDblFinder.class))

# ---- Consensus: flagged by BOTH tools ----
obj@meta.data$consensus_doublet <- (
  obj@meta.data$doublet_finder    == "Doublet" &   # DoubletFinder uses capital D
  obj@meta.data$scDblFinder.class == "doublet"     # scDblFinder uses lowercase d
)
# Label used for the QC violin plot below
obj@meta.data$consensus_doublet_label <- ifelse(obj@meta.data$consensus_doublet %in% TRUE, "Doublet", "Singlet")

# ---- Overlap statistics ----
df_doublet_barcodes   <- rownames(sglt_dblt_metadata)[sglt_dblt_metadata$doublet_finder == "Doublet"]
scdbl_doublet_barcodes <- rownames(meta_scdblfinder)[meta_scdblfinder$scDblFinder.class == "doublet"]
doublet_overlap        <- length(intersect(df_doublet_barcodes, scdbl_doublet_barcodes))

cat("\n=== Doublet detection summary ===\n")
cat("DoubletFinder doublets:        ", length(df_doublet_barcodes), "\n")
cat("scDblFinder doublets:          ", length(scdbl_doublet_barcodes), "\n")
cat("Consensus (overlap) doublets:  ", doublet_overlap, "\n")
cat("Cells retained after filtering:", sum(!obj@meta.data$consensus_doublet, na.rm = TRUE), "\n")

# ---- Per-sample summary ----
doublets_summary <- obj@meta.data %>%
  group_by(orig.ident) %>%
  summarise(
    total           = n(),
    DF_doublets     = sum(doublet_finder    == "Doublet", na.rm = TRUE),
    scDbl_doublets  = sum(scDblFinder.class == "doublet", na.rm = TRUE),
    consensus       = sum(consensus_doublet, na.rm = TRUE),
    pct_removed     = round(100 * sum(consensus_doublet, na.rm = TRUE) / n(), 2),
    .groups = "drop"
  )
print(doublets_summary, n = Inf)
write.table(doublets_summary,
            file  = paste0(plot_dir, "consensus_doublets_summary.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")

# ---- Plot 1: Overlap bar chart ----
doublet_data <- data.frame(
  Method   = c("scDblFinder", "DoubletFinder", "Consensus overlap"),
  Doublets = c(length(scdbl_doublet_barcodes),
               length(df_doublet_barcodes),
               doublet_overlap)
)

p_overlap <- ggplot(doublet_data, aes(x = Method, y = Doublets, fill = Method)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.8) +
  scale_fill_manual(values = c("#4C72B0", "#DD8452", "#55A868")) +
  geom_text(aes(label = Doublets), vjust = -0.4, size = 6, fontface = "bold") +
  labs(title = "Doublets identified by each method",
       y = "Number of doublets", x = "") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"),
        panel.grid.minor   = element_blank()) +
  ylim(0, max(doublet_data$Doublets) * 1.2)

ggsave(paste0(plot_dir, "doublet_methods_overlap_bar.png"),
       p_overlap, dpi = 300, height = 8, width = 10)

# ---- Plot 2: Per-sample stacked bar (DoubletFinder) ----
df_per_sample <- obj@meta.data %>%
  group_by(orig.ident, doublet_finder) %>%
  summarise(total_count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(countT  = sum(total_count),
         percent = paste0(round(100 * total_count / countT, 2), "%"),
         percent_label = ifelse(doublet_finder == "Doublet", percent, NA),
         label_position = ifelse(doublet_finder == "Doublet", countT, NA)) %>%
  ungroup() %>%
  mutate(doublet_finder = factor(doublet_finder, levels = c("Doublet", "Singlet")))

p_df <- ggplot(df_per_sample, aes(x = orig.ident, y = total_count, fill = doublet_finder)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black") +
  geom_text(aes(y = label_position + 200, label = percent_label),
            size = 3, fontface = "bold", na.rm = TRUE) +
  scale_fill_manual(values = c("Doublet" = "#F44336", "Singlet" = "#4CAF50"),
                    name = "") +
  labs(title = "DoubletFinder: Singlets and Doublets per Sample",
       x = "Sample", y = "Cell Count") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(hjust = 0.5, face = "bold"),
        legend.position = "top")

ggsave(paste0(plot_dir, "DoubletFinder_per_sample_bar.png"),
       p_df, dpi = 300, height = 8, width = 16)

# ---- Plot 3: Per-sample stacked bar (scDblFinder) ----
scdbl_per_sample <- obj@meta.data %>%
  group_by(orig.ident, scDblFinder.class) %>%
  summarise(total_count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(countT  = sum(total_count),
         percent = paste0(round(100 * total_count / countT, 2), "%"),
         percent_label  = ifelse(scDblFinder.class == "doublet", percent, NA),
         label_position = ifelse(scDblFinder.class == "doublet", countT, NA)) %>%
  ungroup() %>%
  mutate(scDblFinder.class = factor(scDblFinder.class, levels = c("doublet", "singlet")))

p_scdbl <- ggplot(scdbl_per_sample, aes(x = orig.ident, y = total_count, fill = scDblFinder.class)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black") +
  geom_text(aes(y = label_position + 200, label = percent_label),
            size = 3, fontface = "bold", na.rm = TRUE) +
  scale_fill_manual(values = c("doublet" = "#F44336", "singlet" = "#4CAF50"),
                    name = "") +
  labs(title = "scDblFinder: Singlets and Doublets per Sample",
       x = "Sample", y = "Cell Count") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(hjust = 0.5, face = "bold"),
        legend.position = "top")

ggsave(paste0(plot_dir, "scDblFinder_per_sample_bar.png"),
       p_scdbl, dpi = 300, height = 8, width = 16)

# ---- Plot 4: QC violin comparing singlets vs consensus doublets ----
# Extract data manually
df <- FetchData(obj, vars = c("nFeature_RNA", "percent.mt",
                              "orig.ident", "consensus_doublet_label"),
                layer = "counts")

# Violin plot function
vln_plot <- function(yvar) {
  ggplot(df, aes(x = orig.ident, y = .data[[yvar]], fill = consensus_doublet_label)) +
    geom_violin(scale = "width", position = position_dodge(0.9)) +
    scale_fill_manual(values = c("Doublet" = "#F44336", "Singlet" = "#4CAF50"),
                      name = "") +
    theme_classic() +
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = yvar)
}

p_vln <- vln_plot("nFeature_RNA") / vln_plot("percent.mt")

ggsave(paste0(plot_dir, "consensus_doublets_QC_violin.png"),
       p_vln, dpi = 300, height = 10, width = 20)

# ---- Filter consensus doublets ----
obj <- subset(obj, subset = consensus_doublet == FALSE)
cat("\nFinal cell count after doublet removal:", ncol(obj), "\n")

# ---- Save ----
qs_save(obj, paste0(save_dir, "obj_merge_filt_doublets_removed.qs2"))
cat("Saved: obj_merge_filt_doublets_removed.qs2\n")
