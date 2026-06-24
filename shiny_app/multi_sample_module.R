# Multi-sample integration workflow.

# Design: the user runs the full single-sample pipeline on each sample, then uploads both clustered RDS files here. The worker merges on shared genes, log-normalises, runs Harmony batch correction, 
# builds a joint UMAP + clustering, and runs DESeq2 pseudobulk differential expression between the two conditions.

# Landing card: upload Sample A / Sample B (.rds) with editable display labels.
multiSampleLandingUI <- function() { # [Redacted — UI layout] }

# Compare-tab sidebar: integration settings (Harmony vs merge-only, PCs, resolution), joint-UMAP point size, pseudobulk-DE cluster selector, and exports.
multiSampleControlsUI <- function() { # [Redacted — UI layout] }

# Compare-tab main panel: joint UMAP by sample and by cluster, cluster-composition bar chart, and the differential-expression results section.
multiSampleExploreUI <- function() { # [Redacted — UI layout] }

# Multi-sample server: handles sample uploads/validation; submits the integrate_samples and pseudobulk_de jobs to the worker; renders the joint UMAPs, composition plot, and DE volcano/table 
  # (with a log2FC-only fallback when there is one replicate per condition), and provides integrated-object and plot exports.
multiSampleServer <- function(input, output, session,
                              filtered_obj, update_filtered_obj, pipeline_log,
                              active_job, queue_dir_rv, results_dir_rv,
                              SESSION_ID, show_landing, raw_obj,
                              ms_mode, ms_sample_a, ms_sample_b,
                              ms_integrated, ms_de_results,
                              bin_size = reactive("square_008um")) {
  # [Redacted — upload handling, job submission, integrated plots, DE results]
}
