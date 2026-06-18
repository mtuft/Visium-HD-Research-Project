# 01_cross_sample_overview.R
# Whole-sample comparison of Phase 1 vs BANKSY clustering across P1, P2, P5.
# Produces: cross_sample_summary.csv per-sample summary CSV, cell-type enrichment CSV, spatial comparison PDF, cross_sample_cluster_comparison.pdf

source("/scratch/prj/mmg_grp_single_cell_lab/projects/michael/pipeline_output/CRC_Comparative_Analysis/scripts/00_helpers.R")
library(ggplot2)
library(patchwork)
library(pheatmap)
library(dplyr)
library(tidyr)

ensure_dirs()

analyze_sample <- function(sample_name) {
  message(sample_name)
  sample_dir <- file.path(OUTPUT_DIR, sample_name)
  
  phase1 <- readRDS(phase1_path(sample_name))
  banksy <- readRDS(banksy_path(sample_name))
  
  n_phase1 <- length(unique(phase1$seurat_clusters))
  n_banksy <- length(unique(banksy$seurat_clusters))
  total_spots <- ncol(phase1)
  
  # 1. Confusion matrix
  confusion <- table(
    Phase1 = as.character(phase1$seurat_clusters),
    BANKSY = as.character(banksy$seurat_clusters)
  )
  pdf(file.path(sample_dir, "01_confusion_matrix.pdf"), width = 12, height = 10)
  pheatmap(log10(confusion + 1),
           cluster_rows = FALSE, cluster_cols = FALSE,
           main = paste(sample_name, "- Phase 1 vs BANKSY Cluster Overlap"),
           color = colorRampPalette(c("white", "steelblue", "darkblue"))(50))
  dev.off()
  
  # 2. Concordance
  identical_spots <- sum(as.character(phase1$seurat_clusters) ==
                           as.character(banksy$seurat_clusters))
  concordance_pct <- round(100 * identical_spots / total_spots, 1)
  message(sprintf("  Concordance: %.1f%% (%d/%d)",
                  concordance_pct, identical_spots, total_spots))
  
  # 3. BANKSY-specific populations
  banksy_specific <- detect_banksy_specific(phase1, banksy, overlap_threshold = 0.2)
  banksy$population_type <- ifelse(
    as.character(banksy$seurat_clusters) %in% banksy_specific,
    "banksy_specific", "shared"
  )
  n_specific <- sum(banksy$population_type == "banksy_specific")
  pct_specific <- round(100 * n_specific / total_spots, 2)
  
  # 4. Cell-type enrichment in BANKSY-specific vs shared
  if ("deconv_type1" %in% colnames(banksy@meta.data) && length(banksy_specific) > 0) {
    spec_t <- table(banksy$deconv_type1[banksy$population_type == "banksy_specific"])
    shar_t <- table(banksy$deconv_type1[banksy$population_type == "shared"])
    spec_p <- prop.table(spec_t) * 100
    shar_p <- prop.table(shar_t) * 100
    all_types <- unique(c(names(spec_p), names(shar_p)))
    enrichment_df <- data.frame(
      CellType     = all_types,
      Specific_Pct = as.numeric(spec_p[all_types]),
      Shared_Pct   = as.numeric(shar_p[all_types])
    )
    enrichment_df$Enrichment <- enrichment_df$Specific_Pct / enrichment_df$Shared_Pct
    enrichment_df <- enrichment_df[!is.na(enrichment_df$Enrichment), ]
    enrichment_df <- enrichment_df[order(-enrichment_df$Enrichment), ]
    write.csv(enrichment_df,
              file.path(sample_dir, "02_celltype_enrichment.csv"),
              row.names = FALSE)
  }
  
  # 5. Three-panel spatial comparison
  cc <- get_clean_coords(banksy)
  coords <- cc$coords_mat[cc$valid_rows, ]
  sf <- cc$sf
  plot_df <- data.frame(
    x = coords[, 1] * sf,
    y = -coords[, 2] * sf,
    phase1_cluster = as.character(phase1$seurat_clusters[cc$valid_rows]),
    banksy_cluster = as.character(banksy$seurat_clusters[cc$valid_rows]),
    population     = banksy$population_type[cc$valid_rows]
  )
  
  p1 <- ggplot(plot_df, aes(x = x, y = y, color = phase1_cluster)) +
    geom_point(size = 0.06, alpha = 0.8) +
    coord_fixed() + theme_void() +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold")) +
    labs(title = sprintf("Phase 1 (%d clusters)", n_phase1)) +
    guides(color = guide_legend(override.aes = list(size = 2), nrow = 2))
  
  p2 <- ggplot(plot_df, aes(x = x, y = y, color = banksy_cluster)) +
    geom_point(size = 0.06, alpha = 0.8) +
    coord_fixed() + theme_void() +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold")) +
    labs(title = sprintf("BANKSY (%d clusters)", n_banksy)) +
    guides(color = guide_legend(override.aes = list(size = 2), nrow = 3))
  
  p3 <- ggplot(plot_df, aes(x = x, y = y)) +
    geom_point(data = subset(plot_df, population == "shared"),
               color = "gray85", size = 0.03) +
    geom_point(data = subset(plot_df, population == "banksy_specific"),
               aes(color = banksy_cluster), size = 1.2) +
    coord_fixed() + theme_void() +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold")) +
    labs(title = sprintf("BANKSY-specific (n=%d, %.2f%%)", n_specific, pct_specific)) +
    guides(color = guide_legend(override.aes = list(size = 3), nrow = 2))
  
  combined <- (p1 | p2 | p3) +
    plot_annotation(
      title    = paste(sample_name, "- Spatial Clustering Comparison"),
      subtitle = sprintf("Concordance: %.1f%% | BANKSY-specific clusters: %d",
                         concordance_pct, length(banksy_specific))
    )
  ggsave(file.path(sample_dir, "04_spatial_comparison.pdf"),
         combined, width = 18, height = 7, bg = "white")
  
  # 6. Summary row
  summary_stats <- data.frame(
    Sample                = sample_name,
    Total_Spots           = total_spots,
    Phase1_Clusters       = n_phase1,
    BANKSY_Clusters       = n_banksy,
    Concordance_Pct       = concordance_pct,
    BANKSY_Specific_Clusters = length(banksy_specific),
    BANKSY_Specific_Spots = n_specific,
    BANKSY_Specific_Pct   = pct_specific
  )
  write.csv(summary_stats,
            file.path(sample_dir, "00_summary_stats.csv"),
            row.names = FALSE)
  message("  ✓ ", sample_name, " complete")
  summary_stats
}

# Run for all samples
all_summaries <- lapply(SAMPLES, analyze_sample)
combined_summary <- do.call(rbind, all_summaries)

write.csv(combined_summary,
          file.path(OUTPUT_DIR, "cross_sample_summary.csv"),
          row.names = FALSE)
print(combined_summary)

# Cross-sample cluster count comparison plot
plot_data <- combined_summary %>%
  dplyr::select(Sample, Phase1_Clusters, BANKSY_Clusters) %>%
  tidyr::pivot_longer(cols = c(Phase1_Clusters, BANKSY_Clusters),
                      names_to = "Method", values_to = "Clusters")

p_comparison <- ggplot(plot_data, aes(x = Sample, y = Clusters, fill = Method)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = Clusters),
            position = position_dodge(width = 0.9), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Phase1_Clusters" = "steelblue",
                               "BANKSY_Clusters" = "coral"),
                    labels = c("BANKSY", "Phase 1")) +
  labs(title = "Cluster Count Comparison Across Samples",
       y = "Number of Clusters", x = "Sample") +
  theme_minimal() +
  theme(legend.title = element_blank(),
        plot.title   = element_text(face = "bold", size = 14))

ggsave(file.path(OUTPUT_DIR, "cross_sample_cluster_comparison.pdf"),
       p_comparison, width = 8, height = 6)

message("Outputs in:", OUTPUT_DIR)
