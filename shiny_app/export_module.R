# export_module.R
# Figure export, methods paragraph, object/marker export, pathway analysis, GSEA.

# Explore-tab sidebar: metadata filters and the Figure Export panel (title/caption, auto-generated methods paragraph, scale-bar/methods-page toggles, spatial/UMAP/violin/all PDF exports, 
# processed-object .rds export, parameter log).
exploreExportUI <- function() { # [Redacted — UI layout] }

# Markers-tab sidebar: marker/heatmap CSV+PDF exports, Loupe Browser cluster CSV, and the reproducibility pipeline-settings JSON.
exportControlsUI <- function() { # [Redacted — UI layout] }

# Pipeline log panel UI.
pipelineLogUI <- function() { # [Redacted — UI layout] }

# Marker-column key UI.
markersKeyUI <- function() { # [Redacted — UI layout] }

# Export server: auto-generates an editable methods paragraph from the run parameters; renders the marker heatmap; builds the scale bar and assembles spatial/UMAP/violin PDF exports 
  # (with optional methods page); exports the processed Seurat object, Loupe cluster CSV, marker CSVs, parameter log, and pipeline-settings JSON; and runs clusterProfiler GO/KEGG pathway analysis and
  # fgsea/msigdbr GSEA with their plots, tables, and downloads.
exportServer <- function(input, output, session,
                         filtered_obj, markers_data,
                         current_plot, qc_params_store,
                         pipeline_log, selection_markers,
                         brushed_cells        = reactive(NULL),
                         selected_cluster_val = reactive(NULL),
                         cluster_labels_rv    = reactive(list())) {
  # [Redacted — methods text, figure/object/marker exports, pathway analysis, GSEA]
}
