# 00_helpers.R
# Shared functions used across all CRC Phase 1 vs BANKSY comparison scripts.

library(Seurat)
library(spdep)    
library(jsonlite)  
library(png)   

# File paths
BASE_DIR    <- "/scratch/prj/mmg_grp_single_cell_lab/projects/michael/pipeline_output"
DATA_DIR    <- "/scratch/prj/mmg_grp_single_cell_lab/projects/michael/data/CRC_VisiumHD"
INPUT_DIR   <- file.path(BASE_DIR, "CRC_Comparative_Analysis")        # original analysis files (read-only)
OUTPUT_DIR  <- file.path(BASE_DIR, "CRC_Comparative_Analysis_v2")     # new consolidated outputs
SAMPLES     <- c("P1", "P2", "P5")

phase1_path <- function(sample) {
  file.path(BASE_DIR, "CRC_Phase1",
            paste0("crc_", tolower(sample), "_clustered.rds"))
}

banksy_path <- function(sample) {
  file.path(BASE_DIR, "CRC_Phase2_BANKSY",
            paste0("crc_", tolower(sample), "_banksy_clustered.rds"))
}

# Composition metrics
calculate_purity <- function(prop_matrix) {
  row_maxes <- apply(prop_matrix, 1, max)
  mean(row_maxes[!is.na(row_maxes) & is.finite(row_maxes)])
}

calculate_entropy <- function(prop_matrix) {
  # Mean Shannon entropy per row, computed in proportions (0-1).
  entropies <- apply(prop_matrix / 100, 1, function(row) {
    row <- row[row > 0]
    if (length(row) == 0) return(0)
    -sum(row * log2(row))
  })
  mean(entropies)
}

# Spatial coherence — % of k=8 nearest spatial neighbours in same cluster. Returns a numeric vector of per-spot coherence percentages.
calculate_spatial_coherence <- function(coords_mat, cluster_labels, k = 8) {
  stopifnot(nrow(coords_mat) == length(cluster_labels))
  storage.mode(coords_mat) <- "numeric"
  nb <- knn2nb(knearneigh(coords_mat, k = k))
  sapply(seq_len(nrow(coords_mat)), function(i) {
    neighbours <- nb[[i]]
    same <- sum(cluster_labels[neighbours] == cluster_labels[i])
    100 * same / length(neighbours)
  })
}

# Coordinate extraction with NA filtering. Returns list with coords_mat (numeric matrix of x,y), valid_rows (logical vector), and sf (plotting scale factor).
get_clean_coords <- function(seurat_obj) {
  coords <- GetTissueCoordinates(seurat_obj)
  coords_mat <- if (ncol(coords) == 2) as.matrix(coords) else as.matrix(coords[, 1:2])
  storage.mode(coords_mat) <- "numeric"
  valid_rows <- complete.cases(coords_mat)
  sf <- seurat_obj@images[[1]]@scale.factors$hires / 10
  list(coords_mat = coords_mat, valid_rows = valid_rows, sf = sf)
}

# Identify BANKSY-specific clusters: BANKSY clusters that don't overlap with any single Phase 1 cluster by >threshold (default 20%). 
detect_banksy_specific <- function(phase1_obj, banksy_obj, overlap_threshold = 0.2) {
  confusion <- table(
    Phase1 = as.character(phase1_obj$seurat_clusters),
    BANKSY = as.character(banksy_obj$seurat_clusters)
  )
  banksy_clusters <- colnames(confusion)
  banksy_specific <- banksy_clusters[
    sapply(banksy_clusters, function(bc) {
      max(confusion[, bc]) / sum(confusion[, bc]) < overlap_threshold
    })
  ]
  banksy_specific
}

# Subcluster a Phase 1 cluster and return its composition matrices for both methods. 
# Used by 03 and 04. Returns list with standard_prop, banksy_prop, n_standard, n_banksy_significant, plus the raw subcluster assignments for downstream spatial plotting.

subcluster_and_compose <- function(phase1_obj, banksy_obj, cluster_id,
                                   resolution = 0.8, min_banksy_spots = 100,
                                   verbose = FALSE) {
  cluster_spots <- which(phase1_obj$seurat_clusters == cluster_id)
  cluster_obj <- phase1_obj[, cluster_spots]
  
  # Standard subclustering (fixed seed for reproducibility)
  cluster_obj <- NormalizeData(cluster_obj, verbose = FALSE)
  cluster_obj <- FindVariableFeatures(cluster_obj, nfeatures = 3000, verbose = FALSE)
  cluster_obj <- ScaleData(cluster_obj, verbose = FALSE)
  cluster_obj <- RunPCA(cluster_obj, npcs = 20, verbose = FALSE)
  set.seed(42)
  cluster_obj <- FindNeighbors(cluster_obj, dims = 1:20, verbose = FALSE)
  cluster_obj <- FindClusters(cluster_obj, resolution = resolution, verbose = FALSE)
  
  standard_sub <- cluster_obj$seurat_clusters
  standard_ct  <- cluster_obj$deconv_type1
  standard_comp <- table(standard_sub, standard_ct)
  standard_prop <- prop.table(standard_comp, margin = 1) * 100
  
  # BANKSY subclusters for same spots, filtered to those with >min_banksy_spots
  banksy_sub_all <- banksy_obj$seurat_clusters[cluster_spots]
  banksy_ct_all  <- banksy_obj$deconv_type1[cluster_spots]
  sig_banksy <- names(which(table(banksy_sub_all) > min_banksy_spots))
  keep <- banksy_sub_all %in% sig_banksy
  
  banksy_sub <- factor(as.character(banksy_sub_all[keep]))
  banksy_ct  <- factor(as.character(banksy_ct_all[keep]))
  banksy_comp <- table(banksy_sub, banksy_ct)
  banksy_prop <- prop.table(banksy_comp, margin = 1) * 100
  
  if (verbose) {
    message(sprintf("  Cluster %s: %d spots", cluster_id, length(cluster_spots)))
    message(sprintf("  Standard: %d subclusters", length(unique(standard_sub))))
    message(sprintf("  BANKSY:   %d significant subclusters (>%d spots)",
                    length(sig_banksy), min_banksy_spots))
  }
  
  list(
    cluster_spots         = cluster_spots,
    standard_subclusters  = standard_sub,
    standard_composition  = standard_comp,
    standard_prop         = standard_prop,
    n_standard            = length(unique(standard_sub)),
    banksy_subclusters_all = banksy_sub_all,
    banksy_keep           = keep,
    banksy_subclusters    = banksy_sub,
    banksy_composition    = banksy_comp,
    banksy_prop           = banksy_prop,
    n_banksy_significant  = length(sig_banksy),
    significant_banksy_ids = sig_banksy
  )
}

# Ensure output directories exist
ensure_dirs <- function() {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  for (s in SAMPLES) {
    dir.create(file.path(OUTPUT_DIR, s), recursive = TRUE, showWarnings = FALSE)
  }
}
