# global.R — shared constants, helper functions, and app CSS
# Implementation bodies removed to protect proprietary logic; signatures, constants, and descriptions retained.

# Example dataset path — resolved relative to the app directory at runtime
INPUT_RDS         <- file.path("data", "example_dataset.rds")  # [path parameterised]
DEFAULT_DATA_PATH <- INPUT_RDS
PIPELINE_SEED     <- 42L

# Memory & parallel safety
future::plan("sequential")
options(future.globals.maxSize = 32 * 1024^3)
options(mc.cores = 1)
options(Seurat.object.validate = FALSE)

# Dataset-size thresholds governing the large-data pathways
MAX_CELLS_PIPELINE <- 1000000L
MAX_CELLS_WARNING  <-  800000L
BPCELLS_THRESHOLD  <-  150000L   # convert counts to on-disk above this
SKETCH_THRESHOLD   <-  300000L   # use sketch-based clustering above this
SKETCH_N_CELLS     <-   50000L   # target sketch size

# BPCells availability check (cached at startup; app degrades gracefully if absent).
BPCELLS_AVAILABLE <- tryCatch(requireNamespace("BPCells", quietly = TRUE),
                              error = function(e) FALSE)

# Is this matrix already an on-disk BPCells matrix?
is_bpcells_matrix <- function(mat) { # [Redacted] }

# Convert a counts layer to BPCells on-disk format (session-scoped temp dir; returns the object unchanged if BPCells is unavailable or conversion fails).
convert_to_bpcells <- function(obj, assay, session_id = "default", log_fn = message) {
  # [Redacted — on-disk conversion]
}

# Detect the spatial assay (Spatial > RNA > first available).
detect_spatial_assay <- function(obj) { # [Redacted] }

# Memory diagnostics: estimate object size and warn on large datasets.
estimate_object_size_gb <- function(obj) { # [Redacted] }
check_object_safety <- function(obj, max_cells = 50000) { # [Redacted] }

# Deterministic reset of reductions/graphs/clustering for re-processing.
reset_seurat_state <- function(obj) { # [Redacted] }

# HTML log-entry formatter for the pipeline log panel.
log_entry <- function(msg, level = "info") { # [Redacted] }

# Null-coalescing operator
`%||%` <- function(a, b) if (!is.null(a)) a else b

# App CSS — steel-blue theme (header banner, step badges, section cards, status boxes, pipeline-log panel, buttons, sidebar, responsive QC plots).
APP_CSS <- HTML("/* [Redacted — themed stylesheet] */")

# Detect which predicted.* columns Azimuth wrote to metadata (shared by analysis_module and azimuth_module).
detect_azimuth_cols <- function(meta) { # [Redacted] }
