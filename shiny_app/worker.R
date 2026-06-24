#!/usr/bin/env Rscript
# worker.R — Visium HD Shiny App background compute worker
# Implementation bodies removed to protect proprietary logic; the header, structure, job-handler signatures, and input/output descriptions are retained to document the design. NOT runnable.

# The "compute process" in the shinyqueue pattern: run independently of Shiny (Rscript worker.R <queue_dir> <results_dir>), it watches the queue directory
# for jobs submitted by the app, dispatches each to the matching handler, and writes result RDS files the app polls for. One worker can serve many concurrent
# sessions; multiple workers can share one queue for parallelism.

# Arguments, parallelism, optional-package detection 
# [Redacted — reads queue/results dirs from args or env vars; sets the future plan from SLURM_CPUS_PER_TASK; detects BPCells / Banksy / Azimuth availability]

PIPELINE_SEED <- 42L   # must match global.R

# Helpers; detects the spatial assay (Spatial > RNA > first available).
detect_spatial_assay <- function(obj) { # [REDACTED] }

# Write a status string the Shiny UI polls (progress, "complete", "ERROR: …").
write_status <- function(results_dir, session_id, job_id, msg) { # [REDACTED] }

# Check whether the user cancelled the job via the UI; abort cleanly if so.
check_cancelled <- function(job_file) { # [REDACTED] }

# Job handler: filter_pca 
# In: session_id, input_rds, min/max counts & features, mito_enabled, max_mt.
# Out: {results_dir}/{session_id}/pca.rds
# Loads the object, applies QC filters, optionally converts counts to BPCells on-disk storage, log-normalises, finds variable features, scales, runs PCA,
# and persists QC parameters for the methods paragraph.
processPCA <- function(job_file, outfolder) {
  # [Redacted — filtering, BPCells conversion, normalisation, PCA]
}

# Job handler: umap_cluster 
# In: session_id, pca_rds, n_pcs_use, resolution, clustering_method.
# Out: {results_dir}/{session_id}/clustered.rds
# Two clustering pathways: a default leverage-score sketch route (subsample, cluster, project labels + UMAP back to all spots) and an opt-in spatially-aware
# BANKSY route (neighbourhood-augmented matrix, PCA, SNN + Louvain), with automatic fallback to sketch when BANKSY is unavailable or the dataset exceeds
# a memory-safety spot limit.
processUMAP <- function(job_file, outfolder) {
  # [Redacted — sketch and BANKSY clustering pathways + label projection]
}

# Job handler: azimuth 
# In: session_id, clustered_rds, reference.  Out: {…}/azimuth.rds
# Azimuth reference mapping against the chosen atlas (downloads the reference on first use), using the primary spatial/RNA assay as query.
processAzimuth <- function(job_file, outfolder) {
  # [Redacted — Azimuth label transfer]
}

# Job handler: find_markers 
# In: session_id, clustered_rds, marker_type ("all"|"selection"), selection_cells.
# Out: {…}/markers.rds or selection_markers.rds
# FindAllMarkers across clusters, or a selection-vs-background contrast.
processMarkers <- function(job_file, outfolder) {
  # [Redacted — differential expression]
}

# Job handler: integrate_samples 
# In: session_id, sample_a_rds, sample_b_rds, labels, method, n_pcs, resolution.
# Out: {…}/integrated.rds
# Merges two clustered samples on shared genes, normalises, runs PCA, applies Harmony batch correction (optional), then joint UMAP + clustering.
processIntegration <- function(job_file, outfolder) {
  # [Redacted — merge, integration, joint clustering]
}

# Job handler: pseudobulk_de 
# In: session_id, input_rds, cluster, sample labels.  Out: {…}/pseudobulk_de.rds
# Pseudobulks a cluster per sample and runs DESeq2 when replicates allow, falling back to log2 fold-change ranking when there is one sample per group.
processPseudobulkDE <- function(job_file, outfolder) {
  # [Redacted — pseudobulk aggregation + DESeq2 / log2FC fallback]
}

# Job handler: find_subclusters 
# In: session_id, input_rds, cluster, graph_name, resolution.
# Out: {…}/subclustered.rds.  FindSubCluster (Louvain) on a chosen cluster.
processSubcluster <- function(job_file, outfolder) {
  # [Redacted — subclustering]
}

# Job handler: rctd 
# In: session_id, input_rds, reference_rds, cell_type_col, rctd_mode, max_cores.
# Out: {…}/rctd.rds.  Builds SpaceXr SpatialRNA + Reference objects (per-type
# downsampling, minimum-cell filtering) and runs RCTD deconvolution.
processRCTD <- function(job_file, outfolder) {
  # [Redacted — RCTD deconvolution]
}

# Dispatch table + lurk loop 
# Maps each job type to its handler, then polls the queue every 2s: claims any queued job (status -> running), dispatches to the handler, and marks it complete.
processFunctions <- list(
  filter_pca = processPCA,             umap_cluster      = processUMAP,
  azimuth    = processAzimuth,         find_markers      = processMarkers,
  integrate_samples = processIntegration, find_subclusters = processSubcluster,
  rctd       = processRCTD,            pseudobulk_de     = processPseudobulkDE
)
# [Redacted — polling loop]
