# annotation.R
# Cell type annotation via PanglaoDB marker-database lookup.

# Scores each cluster's top marker genes against PanglaoDB's curated cell-type marker sets using Fisher's exact test, returning ranked candidate cell types.
# The PanglaoDB table is downloaded once and cached locally (panglaodb_cache.rds); subsequent launches use the cache with no internet required.

PANGLAODB_URL   <- "https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz"
PANGLAODB_CACHE <- NULL   # [Redacted — resolves a cache path in the app directory]

# Load (or download + cache) the PanglaoDB marker table as a named list: cell_type -> character vector of marker gene symbols (human-relevant entries).
load_panglaodb <- function(cache_path = PANGLAODB_CACHE) {
  # [Redacted — download / cache / filter]
}

# Score one cluster's marker genes against the database via Fisher's exact test; returns a data.frame of candidate cell types ranked by adjusted p-value.
score_panglaodb <- function(query_genes, db, universe_n = 20000L) {
  # [Redacted — Fisher's exact test scoring + BH adjustment]
}

# UI card for the PanglaoDB lookup panel.
annotationUI <- function() {
  # [Redacted — UI layout]
}

# Server: loads PanglaoDB once per session; computes top markers for the selected cluster; runs the lookup on demand; renders ranked candidate cell types with
# one-click "Apply" to set the cluster label.
annotationServer <- function(input, output, session,
                             markers_data, selected_cluster_val,
                             cluster_labels_rv, filtered_obj) {
  # [Redacted — lookup orchestration + results UI]
}
