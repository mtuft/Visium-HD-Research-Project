#!/usr/bin/env Rscript
# CRC Phase 2: BANKSY Spatial Clustering

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop("Usage: Rscript banksy_cluster.R <data_dir> <deconv_file> <sample_id> <output_rds> <plot_dir>")
}

DATA_DIR <- args[1]
DECONV_FILE <- args[2]
SAMPLE_ID <- args[3]
OUTPUT_RDS <- args[4]
PLOT_DIR <- args[5]

cat("CRC BANKSY Pipeline -", SAMPLE_ID, "\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Libraries
library(Seurat)
library(SeuratObject)
library(SeuratWrappers)
library(Banksy)
library(ggplot2)
library(dplyr)

# LOAD DATA
cat("[1/9] Loading data...\n")
obj <- Load10X_Spatial(
  data.dir = DATA_DIR,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE
)
cat("Loaded:", ncol(obj), "spots,", nrow(obj), "genes\n\n")

# QC FILTERING
cat("[2/9] Setting QC thresholds...\n")
MIN_COUNTS <- 100
MIN_FEATURES <- 50

cat(sprintf("Minimum nCount:   %d\n", MIN_COUNTS))
cat(sprintf("Minimum nFeature: %d\n", MIN_FEATURES))

cat("[3/9] Filtering spots...\n")
meta <- obj@meta.data
count_col <- "nCount_Spatial"
feat_col <- "nFeature_Spatial"

keep <- (meta[[count_col]] >= MIN_COUNTS &
         meta[[feat_col]] >= MIN_FEATURES)

n_before <- ncol(obj)
n_after <- sum(keep)
cat(sprintf("Before: %d spots\n", n_before))
cat(sprintf("After:  %d spots (%.1f%% retained)\n\n", n_after, 100 * n_after / n_before))

filt <- subset(obj, cells = rownames(meta)[keep])
rm(obj); gc()

# NORMALIZATION
cat("[4/9] Normalizing data...\n")
filt <- NormalizeData(filt, normalization.method = "LogNormalize", scale.factor = 10000)

# VARIABLE FEATURES
cat("[5/9] Finding variable features...\n")
filt <- FindVariableFeatures(filt, selection.method = "vst", nfeatures = 3000)

# SCALING
cat("[6/9] Scaling data...\n")
filt <- ScaleData(filt)

# PCA
cat("[7/9] Running PCA...\n")
filt <- RunPCA(filt, features = VariableFeatures(object = filt), npcs = 20)
cat("PCA complete\n\n")

# BANKSY
cat("[8/9] Running BANKSY (this will take ~10-15 minutes)...\n\n")

set.seed(42)
filt <- RunBanksy(
  filt,
  lambda = 0.2,
  assay = "Spatial",
  slot = "data",
  features = "variable",
  k_geom = 50,
  spatial_mode = "kNN_median"
)

cat("BANKSY matrix created. Running PCA on BANKSY features...\n")
filt <- RunPCA(filt, assay = "BANKSY", features = rownames(filt[["BANKSY"]]), npcs = 20)

cat("Clustering...\n")
filt <- FindNeighbors(filt, reduction = "pca", dims = 1:20, assay = "BANKSY")
filt <- FindClusters(filt, resolution = 0.8, algorithm = 1)

cat("Running UMAP...\n")
filt <- RunUMAP(filt, reduction = "pca", dims = 1:20, assay = "BANKSY")

n_clusters <- length(unique(filt$seurat_clusters))
cat("BANKSY clustering complete:", n_clusters, "clusters\n")
cat("Cluster levels:", paste(levels(filt$seurat_clusters), collapse = ", "), "\n\n")

# DECONVOLUTION INTEGRATION
cat("Loading deconvolution data...\n")
deconv_df <- read.csv(DECONV_FILE)
cat("Loaded deconvolution data:", nrow(deconv_df), "entries\n\n")

# Match barcodes
barcode_match <- match(colnames(filt), deconv_df$barcode)

# Add deconvolution metadata
filt$deconv_class <- NA
filt$deconv_class[!is.na(barcode_match)] <- deconv_df$DeconvolutionClass[barcode_match[!is.na(barcode_match)]]

filt$deconv_type1 <- NA
filt$deconv_type1[!is.na(barcode_match)] <- deconv_df$DeconvolutionLabel1[barcode_match[!is.na(barcode_match)]]

filt$deconv_type2 <- NA
filt$deconv_type2[!is.na(barcode_match)] <- deconv_df$DeconvolutionLabel2[barcode_match[!is.na(barcode_match)]]

# Summary
spots_with_deconv <- sum(!is.na(filt$deconv_class))
cat(sprintf("Deconvolution summary (%d spots with predictions):\n", spots_with_deconv))

deconv_summary <- table(filt$deconv_class, useNA = "ifany")
for (class in names(deconv_summary)) {
  cat(sprintf("  %-18s %d \n", paste0(class, ":"), deconv_summary[class]))
}

cat("\nTop 10 cell types detected:\n\n")
top_types <- sort(table(filt$deconv_type1), decreasing = TRUE)[1:10]
print(top_types)
cat("\n")

# SAVE OUTPUTS
cat("[9/9] Saving plots and output...\n")

# IGV color palette
igv_colors <- c("#5773CC", "#FFB900", "#2FA236", "#CC66FF", "#F97A6D", 
                "#33CCCC", "#FF6600", "#009999", "#993366", "#FFCC00",
                "#996600", "#66CCFF", "#FF9999", "#669933", "#CC99FF",
                "#FF99CC", "#99CCFF", "#CCFF99", "#FF6699", "#FFCC99",
                "#99CCCC", "#CC9999", "#CCCCFF", "#FFCCCC", "#99FF99")

cluster_colors <- igv_colors[1:n_clusters]
names(cluster_colors) <- levels(filt$seurat_clusters)

# UMAP plot
p1 <- DimPlot(filt, reduction = "umap", group.by = "seurat_clusters", 
              cols = cluster_colors, pt.size = 0.5) +
  ggtitle(paste("UMAP - BANKSY Clusters -", SAMPLE_ID)) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        legend.key.size = unit(4, "mm"))

ggsave(file.path(PLOT_DIR, "04_umap_banksy.pdf"), p1, width = 10, height = 8)

# Spatial plot
p2 <- SpatialDimPlot(filt, group.by = "seurat_clusters", 
                     cols = cluster_colors, pt.size.factor = 1.2) +
  ggtitle(paste("Spatial - BANKSY Clusters -", SAMPLE_ID)) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        legend.key.size = unit(4, "mm"))

ggsave(file.path(PLOT_DIR, "05_spatial_banksy.pdf"), p2, width = 12, height = 10)
ggsave(file.path(PLOT_DIR, "05_spatial_banksy.png"), p2, width = 12, height = 10, dpi = 300)

# Remove BANKSY assay before saving
cat("Removing BANKSY assay to reduce file size...\n")
DefaultAssay(filt) <- "Spatial"
filt[["BANKSY"]] <- NULL

# Save RDS
saveRDS(filt, OUTPUT_RDS)

cat("  -", OUTPUT_RDS, "\n")
cat("  -", file.path(PLOT_DIR, "04_umap_banksy.pdf"), "\n")
cat("  -", file.path(PLOT_DIR, "05_spatial_banksy.pdf"), "\n")
cat("  -", file.path(PLOT_DIR, "05_spatial_banksy.png"), "\n")

# Record session info for reproducibility
writeLines(capture.output(sessionInfo()),
           file.path(PLOT_DIR, "session_info_banksy.txt"))
cat("Clusters:", n_clusters, "\n")
