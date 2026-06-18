# 03_target_cluster.R
# Analysis of three flagship Phase 1 clusters representing distinct spatial architectures:
# P1 cluster 7  (Organised), P2 cluster 6  (Homogeneous), P5 cluster 2  (Intermixed)

# For each cluster, compute:
#   - Standard vs BANKSY subcluster composition (cell-type heatmaps)
#   - Side-by-side spatial layout of subclusters
#   - Per-spot purity mapped back to tissue location (with/without H&E)
#   - kNN spatial coherence at k=10 and k=20

# Outputs (per cluster, in {OUTPUT_DIR}/{sample}/): cluster{N}_composition_standard.pdf, cluster{N}_composition_banksy.pdf, cluster{N}_subcluster_spatial.png,
# cluster{N}_purity_spatial.png, cluster{N}_purity_spatial_with_histology.png, cluster{N}_coherence_violin.png

# Top-level outputs in {OUTPUT_DIR}/:
#   target_cluster_metrics.csv             (long format, one row per method)
#   target_cluster_coherence_per_spot.csv  (long format, per-spot values)

# Coherence definition: for each spot, the fraction of its k nearest spatial neighbours (excluding self) that share its subcluster label.
# Computed via FNN::get.knn on tissue x/y coordinates.

source("/scratch/prj/mmg_grp_single_cell_lab/projects/michael/pipeline_output/CRC_Comparative_Analysis/scripts/00_helpers.R")
library(ggplot2)
library(patchwork)
library(pheatmap)
library(dplyr)
library(FNN)
library(png)
library(jsonlite)
ensure_dirs()

# Configuration
TARGET_CLUSTERS <- tibble::tribble(
  ~Sample, ~Cluster, ~Pattern,
  "P1",    "7",      "Organised",
  "P2",    "6",      "Homogeneous",
  "P5",    "2",      "Intermixed"
)

COHERENCE_K_VALUES  <- c(10, 20)
SHOW_HISTOLOGY      <- TRUE
COMPOSITION_MIN_PCT <- 5  # cell types kept if max % across subclusters > this
RAW_DATA_BASE       <- "/scratch/prj/mmg_grp_single_cell_lab/projects/michael/data/CRC_VisiumHD"

# Helper: kNN spatial coherence
# For each spot, find k nearest spatial neighbours (excluding self), then return the fraction of those neighbours sharing its subcluster label.
compute_knn_coherence <- function(coords, labels, k = 10) {
  stopifnot(nrow(coords) == length(labels))
  if (nrow(coords) <= k + 1) {
    warning(sprintf("Only %d spots; coherence at k=%d undefined", nrow(coords), k))
    return(rep(NA_real_, nrow(coords)))
  }
  labels_chr <- as.character(labels)
  nn_idx     <- FNN::get.knn(coords, k = k)$nn.index
  vapply(seq_len(nrow(coords)), function(i) {
    mean(labels_chr[nn_idx[i, ]] == labels_chr[i])
  }, numeric(1))
}

# Robust max that returns NA for empty/all-NA rows instead of -Inf
robust_max <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else max(x)
}

# Helper: load H&E histology image and scale factor for a sample. Returns NULL if either file is missing (caller falls back to clean plot).

load_histology <- function(sample_name) {
  json_path <- file.path(RAW_DATA_BASE, sample_name,
                         "binned_outputs/square_008um/spatial/scalefactors_json.json")
  img_path  <- file.path(RAW_DATA_BASE, sample_name,
                         "binned_outputs/square_008um/spatial/tissue_hires_image.png")
  if (!file.exists(json_path) || !file.exists(img_path)) {
    message("  [histology] Missing files for ", sample_name, " - skipping overlay")
    return(NULL)
  }
  list(img    = png::readPNG(img_path),
       scalef = jsonlite::fromJSON(json_path)$tissue_hires_scalef)
}

# Helper: cell-type composition heatmap
plot_composition_heatmap <- function(prop_matrix, title_text, output_path,
                                     min_pct = COMPOSITION_MIN_PCT) {
  abundant <- names(which(apply(prop_matrix, 2, max, na.rm = TRUE) > min_pct))
  if (length(abundant) == 0) abundant <- colnames(prop_matrix)
  filtered <- prop_matrix[, abundant, drop = FALSE]
  filtered[is.na(filtered)] <- 0
  
  pdf(output_path, width = 14, height = 10)
  pheatmap(filtered,
           cluster_rows = TRUE, cluster_cols = TRUE,
           main         = title_text,
           color        = colorRampPalette(c("white", "orange", "red", "darkred"))(50),
           display_numbers = TRUE, number_format = "%.1f",
           fontsize_number = 7, fontsize_row = 10, fontsize_col = 9,
           angle_col    = 45,
           cellwidth    = 20, cellheight = 20,
           border_color = "grey80")
  dev.off()
}

# Helper: side-by-side spatial layout of subclusters
plot_subcluster_spatial <- function(coords_std, labels_std,
                                    coords_bks, labels_bks,
                                    n_std, n_bks_sig,
                                    sample_name, cluster_id,
                                    output_path, sf) {
  df_std <- data.frame(x = coords_std[, 1] * sf, y = -coords_std[, 2] * sf,
                       Subcluster = factor(labels_std))
  df_bks <- data.frame(x = coords_bks[, 1] * sf, y = -coords_bks[, 2] * sf,
                       Subcluster = factor(labels_bks))
  
  make_panel <- function(df, title_text) {
    ggplot(df, aes(x, y, color = Subcluster)) +
      geom_point(size = 1.0, alpha = 0.8) +
      coord_fixed() + theme_void() +
      labs(title = title_text) +
      theme(plot.title      = element_text(size = 12, face = "bold", hjust = 0.5),
            legend.position = "right",
            legend.key.size = unit(0.4, "cm"),
            legend.text     = element_text(size = 8))
  }
  
  p_std <- make_panel(df_std,
                      sprintf("Standard Subclustering\n%d subclusters", n_std))
  p_bks <- make_panel(df_bks,
                      sprintf("BANKSY Subclustering\n%d subclusters (>100 spots)",
                              n_bks_sig))
  
  combined <- p_std + p_bks +
    plot_annotation(
      title = sprintf("%s Cluster %s: Subcluster Spatial Organisation",
                      sample_name, cluster_id),
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  ggsave(output_path, combined, width = 14, height = 7, dpi = 150)
}

# Helper: per-spot purity heatmap (clean background or H&E overlay)
plot_purity_spatial <- function(coords_std, purity_std,
                                coords_bks, purity_bks,
                                n_std, n_bks_sig,
                                mean_purity_std, mean_purity_bks,
                                sample_name, cluster_id,
                                output_path, sf, hist = NULL) {
  
  if (!is.null(hist)) {
    sc <- hist$scalef
    df_std <- data.frame(x = coords_std[, 1] * sc, y = -coords_std[, 2] * sc,
                         purity = purity_std)
    df_bks <- data.frame(x = coords_bks[, 1] * sc, y = -coords_bks[, 2] * sc,
                         purity = purity_bks)
    all_x <- c(df_std$x, df_bks$x); all_y <- c(df_std$y, df_bks$y)
    pad_x <- 0.05 * (max(all_x) - min(all_x))
    pad_y <- 0.05 * (max(all_y) - min(all_y))
    x_min <- max(1, floor(min(all_x) - pad_x))
    x_max <- min(ncol(hist$img), ceiling(max(all_x) + pad_x))
    y_min <- max(1, floor(min(-all_y) - pad_y))
    y_max <- min(nrow(hist$img), ceiling(max(-all_y) + pad_y))
    img_crop <- hist$img[y_min:y_max, x_min:x_max, ]
    raster_layer <- annotation_raster(img_crop,
                                      xmin = x_min, xmax = x_max,
                                      ymin = -y_max, ymax = -y_min)
  } else {
    df_std <- data.frame(x = coords_std[, 1] * sf, y = -coords_std[, 2] * sf,
                         purity = purity_std)
    df_bks <- data.frame(x = coords_bks[, 1] * sf, y = -coords_bks[, 2] * sf,
                         purity = purity_bks)
    raster_layer <- NULL
  }
  
  purity_lim <- range(c(df_std$purity, df_bks$purity), na.rm = TRUE)
  
  make_panel <- function(df, title_text) {
    p <- ggplot()
    if (!is.null(raster_layer)) p <- p + raster_layer
    p +
      geom_point(data = df, aes(x, y, color = purity),
                 size = 1.4, alpha = 0.85) +
      scale_color_viridis_c(name   = "Cell-type\nPurity (%)",
                            limits = purity_lim,
                            option = "plasma",
                            breaks = seq(20, 100, 20)) +
      coord_fixed() + theme_void() +
      labs(title = title_text) +
      theme(plot.title      = element_text(size = 12, face = "bold", hjust = 0.5),
            legend.position = "right",
            legend.title    = element_text(size = 10),
            legend.text     = element_text(size = 8))
  }
  
  p_std <- make_panel(df_std,
                      sprintf("Standard Subclustering\n%d subclusters, %.1f%% purity",
                              n_std, mean_purity_std))
  p_bks <- make_panel(df_bks,
                      sprintf("BANKSY Subclustering\n%d subclusters, %.1f%% purity",
                              n_bks_sig, mean_purity_bks))
  
  combined <- p_std + p_bks +
    plot_annotation(
      title = sprintf("%s Cluster %s: Spatial Cell-Type Purity",
                      sample_name, cluster_id),
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  ggsave(output_path, combined, width = 14, height = 7, dpi = 150)
}

# Helper: per-spot coherence violin plot, faceted by k
plot_coherence_violin <- function(coh_std_k10, coh_bks_k10,
                                  coh_std_k20, coh_bks_k20,
                                  sample_name, cluster_id, output_path) {
  df <- bind_rows(
    data.frame(Method = "Standard", k = "k = 10", coherence = coh_std_k10 * 100),
    data.frame(Method = "BANKSY",   k = "k = 10", coherence = coh_bks_k10 * 100),
    data.frame(Method = "Standard", k = "k = 20", coherence = coh_std_k20 * 100),
    data.frame(Method = "BANKSY",   k = "k = 20", coherence = coh_bks_k20 * 100)
  )
  df$Method <- factor(df$Method, levels = c("Standard", "BANKSY"))
  
  p <- ggplot(df, aes(x = Method, y = coherence, fill = Method)) +
    geom_violin(alpha = 0.5, scale = "width") +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white",
                 colour = "#333") +
    scale_fill_manual(values = c("Standard" = "#4a7fa5",
                                 "BANKSY"   = "#c0392b")) +
    facet_wrap(~ k) +
    labs(title    = sprintf("%s Cluster %s: kNN Spatial Coherence",
                            sample_name, cluster_id),
         subtitle = "% of nearest spatial neighbours sharing subcluster label",
         x = NULL, y = "Coherence (%)") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none",
          plot.title      = element_text(face = "bold", colour = "#2c5f7a"),
          plot.subtitle   = element_text(colour = "#666", size = 10))
  ggsave(output_path, p, width = 9, height = 6, dpi = 150)
}

# Main loop
metrics_rows   <- list()
coherence_rows <- list()

for (i in seq_len(nrow(TARGET_CLUSTERS))) {
  row         <- TARGET_CLUSTERS[i, ]
  sample_name <- row$Sample
  cluster_id  <- row$Cluster
  pattern     <- row$Pattern
  
 
  message(sprintf("Cluster",
                  sample_name, cluster_id, pattern))
  sample_dir <- file.path(OUTPUT_DIR, sample_name)
  if (!dir.exists(sample_dir)) dir.create(sample_dir, recursive = TRUE)
  phase1 <- readRDS(phase1_path(sample_name))
  banksy <- readRDS(banksy_path(sample_name))
  
  res <- subcluster_and_compose(phase1, banksy, cluster_id,
                                resolution = 0.8, min_banksy_spots = 100,
                                verbose = TRUE)
  
  # Spatial coords + scale factor
  sf <- phase1@images[[1]]@scale.factors$hires / 10
  all_coords <- GetTissueCoordinates(phase1)
  if (ncol(all_coords) > 2) all_coords <- all_coords[, 1:2]
  all_coords <- as.matrix(all_coords)
  storage.mode(all_coords) <- "numeric"
  
  cluster_spot_ids <- res$cluster_spots
  coords_cluster   <- all_coords[cluster_spot_ids, , drop = FALSE]
  
  # Standard pathway: all valid spots in cluster
  valid_std  <- complete.cases(coords_cluster)
  coords_std <- coords_cluster[valid_std, , drop = FALSE]
  labels_std <- as.character(res$standard_subclusters)[valid_std]
  
  # BANKSY pathway: only spots in subclusters with >100 spots
  banksy_all <- as.character(res$banksy_subclusters_all)
  banksy_counts <- table(banksy_all)
  significant   <- names(banksy_counts[banksy_counts > 100])
  keep_bks      <- banksy_all %in% significant
  coords_bks    <- coords_cluster[keep_bks, , drop = FALSE]
  labels_bks    <- banksy_all[keep_bks]
  valid_bks     <- complete.cases(coords_bks)
  coords_bks    <- coords_bks[valid_bks, , drop = FALSE]
  labels_bks    <- labels_bks[valid_bks]
  
  # Per-subcluster purity
  std_subcluster_purity <- apply(res$standard_prop, 1, robust_max)
  bks_subcluster_purity <- apply(res$banksy_prop,   1, robust_max)
  purity_std_per_spot <- std_subcluster_purity[labels_std]
  purity_bks_per_spot <- bks_subcluster_purity[labels_bks]
  mean_purity_std     <- mean(std_subcluster_purity, na.rm = TRUE)
  mean_purity_bks     <- mean(bks_subcluster_purity, na.rm = TRUE)
  
  # 1. Composition heatmaps
  message("Generating composition heatmaps...")
  plot_composition_heatmap(
    res$standard_prop,
    title_text = sprintf("%s Standard Subclustering: Cluster %s\n(%d subclusters, %.1f%% purity)",
                         sample_name, cluster_id,
                         res$n_standard, mean_purity_std),
    output_path = file.path(sample_dir,
                            sprintf("cluster%s_composition_standard.pdf", cluster_id))
  )
  plot_composition_heatmap(
    res$banksy_prop,
    title_text = sprintf("%s BANKSY Subclustering: Cluster %s\n(%d subclusters, %.1f%% purity)",
                         sample_name, cluster_id,
                         res$n_banksy_significant, mean_purity_bks),
    output_path = file.path(sample_dir,
                            sprintf("cluster%s_composition_banksy.pdf", cluster_id))
  )
  
  # 2. Subcluster spatial layout
  plot_subcluster_spatial(
    coords_std, labels_std, coords_bks, labels_bks,
    n_std = res$n_standard, n_bks_sig = res$n_banksy_significant,
    sample_name = sample_name, cluster_id = cluster_id,
    output_path = file.path(sample_dir,
                            sprintf("cluster%s_subcluster_spatial.png", cluster_id)),
    sf = sf
  )
  
  # 3. Purity spatial heatmap (clean) 
  plot_purity_spatial(
    coords_std, purity_std_per_spot,
    coords_bks, purity_bks_per_spot,
    n_std = res$n_standard, n_bks_sig = res$n_banksy_significant,
    mean_purity_std = mean_purity_std, mean_purity_bks = mean_purity_bks,
    sample_name = sample_name, cluster_id = cluster_id,
    output_path = file.path(sample_dir,
                            sprintf("cluster%s_purity_spatial.png", cluster_id)),
    sf = sf, hist = NULL
  )
  
  # 4. Purity spatial heatmap (with H&E) 
  if (SHOW_HISTOLOGY) {
    hist_data <- load_histology(sample_name)
    if (!is.null(hist_data)) {
      plot_purity_spatial(
        coords_std, purity_std_per_spot,
        coords_bks, purity_bks_per_spot,
        n_std = res$n_standard, n_bks_sig = res$n_banksy_significant,
        mean_purity_std = mean_purity_std, mean_purity_bks = mean_purity_bks,
        sample_name = sample_name, cluster_id = cluster_id,
        output_path = file.path(sample_dir,
                                sprintf("cluster%s_purity_spatial_with_histology.png",
                                        cluster_id)),
        sf = sf, hist = hist_data
      )
    }
  }
  
  # 5. kNN spatial coherence 
  coh_std_k10 <- compute_knn_coherence(coords_std, labels_std, k = 10)
  coh_bks_k10 <- compute_knn_coherence(coords_bks, labels_bks, k = 10)
  coh_std_k20 <- compute_knn_coherence(coords_std, labels_std, k = 20)
  coh_bks_k20 <- compute_knn_coherence(coords_bks, labels_bks, k = 20)
  
  message(sprintf("  Standard k=10: mean %.1f%% (sd %.1f%%)",
                  100 * mean(coh_std_k10, na.rm = TRUE),
                  100 * sd  (coh_std_k10, na.rm = TRUE)))
  message(sprintf("  BANKSY   k=10: mean %.1f%% (sd %.1f%%)",
                  100 * mean(coh_bks_k10, na.rm = TRUE),
                  100 * sd  (coh_bks_k10, na.rm = TRUE)))
  
  plot_coherence_violin(
    coh_std_k10, coh_bks_k10, coh_std_k20, coh_bks_k20,
    sample_name = sample_name, cluster_id = cluster_id,
    output_path = file.path(sample_dir,
                            sprintf("cluster%s_coherence_violin.png", cluster_id))
  )
  
  # 6. Append summary rows 
  metrics_rows[[paste0(sample_name, "_", cluster_id, "_std")]] <- data.frame(
    Sample = sample_name, Cluster = cluster_id, Pattern = pattern,
    Method = "Standard", N_Spots = length(labels_std),
    N_Subclusters = res$n_standard,
    Mean_Purity   = round(mean_purity_std, 2),
    Mean_Entropy  = round(calculate_entropy(res$standard_prop), 3),
    Coherence_k10_Mean = round(100 * mean(coh_std_k10, na.rm = TRUE), 2),
    Coherence_k10_SD   = round(100 * sd  (coh_std_k10, na.rm = TRUE), 2),
    Coherence_k20_Mean = round(100 * mean(coh_std_k20, na.rm = TRUE), 2),
    Coherence_k20_SD   = round(100 * sd  (coh_std_k20, na.rm = TRUE), 2),
    stringsAsFactors = FALSE
  )
  metrics_rows[[paste0(sample_name, "_", cluster_id, "_bks")]] <- data.frame(
    Sample = sample_name, Cluster = cluster_id, Pattern = pattern,
    Method = "BANKSY", N_Spots = length(labels_bks),
    N_Subclusters = res$n_banksy_significant,
    Mean_Purity   = round(mean_purity_bks, 2),
    Mean_Entropy  = round(calculate_entropy(res$banksy_prop), 3),
    Coherence_k10_Mean = round(100 * mean(coh_bks_k10, na.rm = TRUE), 2),
    Coherence_k10_SD   = round(100 * sd  (coh_bks_k10, na.rm = TRUE), 2),
    Coherence_k20_Mean = round(100 * mean(coh_bks_k20, na.rm = TRUE), 2),
    Coherence_k20_SD   = round(100 * sd  (coh_bks_k20, na.rm = TRUE), 2),
    stringsAsFactors = FALSE
  )
  
  coherence_rows[[paste0(sample_name, "_", cluster_id, "_std")]] <- data.frame(
    Sample = sample_name, Cluster = cluster_id, Pattern = pattern,
    Method = "Standard",
    Coherence_k10 = coh_std_k10, Coherence_k20 = coh_std_k20,
    stringsAsFactors = FALSE
  )
  coherence_rows[[paste0(sample_name, "_", cluster_id, "_bks")]] <- data.frame(
    Sample = sample_name, Cluster = cluster_id, Pattern = pattern,
    Method = "BANKSY",
    Coherence_k10 = coh_bks_k10, Coherence_k20 = coh_bks_k20,
    stringsAsFactors = FALSE
  )
  
  # Write running summaries after each cluster
  write.csv(do.call(rbind, metrics_rows),
            file.path(OUTPUT_DIR, "target_cluster_metrics.csv"),
            row.names = FALSE)
  write.csv(do.call(rbind, coherence_rows),
            file.path(OUTPUT_DIR, "target_cluster_coherence_per_spot.csv"),
            row.names = FALSE)
  
  rm(phase1, banksy); gc()
}

message("Outputs in:", OUTPUT_DIR)
