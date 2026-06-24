# Visium HD Spatial Clustering: Comparative Pipeline & Reproducibility

This repository accompanies an MSc research project on 10x Genomics Visium HD
spatial transcriptomics. It contains (1) a reproducible Nextflow pipeline
comparing standard (expression-only) and spatially-aware (BANKSY) clustering on
three colorectal cancer (CRC) samples, (2) a reproducibility-validation pipeline,
and (3) redacted source for an interactive R Shiny analysis application, included
as evidence of the implementation only.

---

## Scope: what is and isn't in this repository

**Fully reproducible:**
- The Nextflow comparative-clustering pipeline (`main.nf`, `nextflow.config`, `submit.sh`)
- The per-sample clustering scripts for both arms (`scripts/phase1_cluster.R`, `scripts/banksy_cluster.R`)
- The cross-sample comparison and deep-dive analysis scripts (`scripts/00`–`03`)
- The reproducibility-validation pipeline (`visiumhd_reproducibility_pipeline.Rmd`)

**Redacted — evidence only, NOT runnable:**
- The interactive Shiny application (`app_redacted/`). The application is proprietary;
  the source has been redacted so that the architecture, breadth, and design are
  documented for assessment while the implementation is not reproducible. Function
  bodies are removed; signatures, structure, and descriptive comments are retained.

**Not included:**
- The raw Visium HD data and deconvolution results — these are third-party and
  available from the original publication (see *Data Availability*).

---

## Repository structure
```
├── README.md
├── LICENSE
├── main.nf                                   # Nextflow pipeline: runs both arms per sample
├── nextflow.config                           # Profiles (SLURM/local), container, resources
├── submit.sh                                 # SLURM submission wrapper
├── visiumhd_reproducibility_pipeline.Rmd     # Reproducibility validation

├── scripts/
│   ├── phase1_cluster.R                      # Standard (expression-only) clustering arm
│   ├── banksy_cluster.R                      # Spatially-aware (BANKSY) clustering arm
│   ├── 00_helpers.R                          # Shared helper functions
│   ├── 01_cross_sample_overview.R            # Cross-sample concordance / coherence
│   ├── 02_highlight_cluster_analysis.R       # (verify name + purpose)
│   └── 03_target_cluster.R                   # (verify name + purpose)

├── shiny_app/
│   ├── global.R                              # Shared constants, helper functions, CSS
│   ├── ui_server.R                           # App UI layout and main server function
│   ├── worker.R                              # Background compute worker (run independently)
│   ├── modules/
│   │   ├── qc_module.R                       # QC filtering, normalisation, PCA
│   │   ├── analysis_module.R                 # Spatial/UMAP plots, marker analysis, volcano plots
│   │   ├── annotation_module.R               # PanglaoDB cell type lookup
│   │   ├── azimuth_module.R                  # Azimuth cell type prediction
│   │   ├── export_module.R                   # Figure export, methods paragraph, save object
│   │   ├── queue_module.R                    # Job queue polling and status management
│   │   └── multi_sample_module.R             # Multi-sample loading and comparison
│   ├── help_tab.R                            # Provides information for users
│   └── app.R                                 # Sources all scripts to run the app

└── env/                                      # sessionInfo() records (add to back §ref in report)

```

---

## App Features

- **Guided QC pipeline** — interactive violin plots with live threshold feedback and real-time spot count
- **Two normalisation methods** — NormalizeData (recommended for spatial) or SCTransform v2
- **PCA with automated PC suggestion** — noise-floor algorithm identifies the elbow; elbow plot with cumulative variance display
- **Standard and BANKSY clustering** — transcriptomics-only (Louvain) or spatially-aware clustering (BANKSY) that incorporates neighbourhood expression context
- **Sketch-based clustering** — leverage-score subsampling for large datasets (>300K spots); auto-enabled above threshold
- **BPCells on-disk storage** — transparent memory management for large datasets (>150K spots); no data is discarded
- **Interactive spatial and UMAP plots** — linked plots, click-to-select clusters, gene expression overlay, toggleable labels
- **Cluster marker analysis** — FindAllMarkers with interactive volcano plots and click-to-explore gene expression
- **Spatial region markers** — draw any spatial selection and run differential expression against the rest of the tissue
- **Pathway analysis** — offline GO and KEGG enrichment via clusterProfiler; no internet required after setup
- **Cell type annotation** — Azimuth reference mapping (13 atlases) and PanglaoDB marker enrichment scoring
- **Custom cluster labelling** — rename clusters with free-text labels; labels propagate across all plots
- **Multi-sample comparison** — side-by-side spatial plots, integrated UMAP, between-sample DE, and cluster composition; batch correction via Harmony or RPCA
- **Publication-ready export** — spatial/UMAP PDFs with scale bars and auto-generated figure captions; methods paragraph and pipeline settings JSON for reproducibility; Loupe Browser CSV export
- **Background job queue** — all heavy computation runs in a separate worker process; the UI remains fully responsive during analysis

---

## The analysis

Each sample (P1, P2, P5) is processed independently through two clustering arms
with identical preprocessing, so any difference in the resulting clusters is
attributable to the spatial augmentation rather than to differing settings:

- **Standard arm** (`phase1_cluster.R`): expression-only clustering. Spots are
  filtered (≥100 UMIs, ≥50 features), log-normalised (scale factor 10,000), 3,000
  variable features, 20 principal components, shared-nearest-neighbour graph,
  Louvain clustering at resolution 0.8.
- **BANKSY arm** (`banksy_cluster.R`): identical preprocessing, with the expression
  matrix augmented by neighbourhood features (BANKSY; λ = 0.2, k_geom = 50,
  `kNN_median`) before PCA and the same SNN + Louvain clustering.

The QC thresholds (≥100 UMIs, ≥50 features) are matched to the source paper's
deconvolution threshold for comparability. Cell-type composition is annotated from
the published RCTD deconvolution results released with the source dataset.

---

## Requirements

Analyses were run in **R 4.4.1** inside the standard **Bioconductor 3.19
Singularity image** (Ubuntu 22.04.4 LTS). Key packages: Seurat v5, BANKSY,
SeuratWrappers, spacexr (RCTD), clusterProfiler. The complete `sessionInfo()` for
each run is in `env/`. Nextflow ≥ 23.10.0 is required. R ≥ 4.2.0 is recommended. The app uses Seurat v5 and relies on the layer-based assay structure introduced in that version.

---

### Required packages

```r
install.packages(c(
  "shiny", "shinyjs", "shinycssloaders",
  "Seurat", "SeuratObject",
  "ggplot2", "ggrepel", "patchwork",
  "dplyr", "stringr",
  "DT", "plotly",
  "future", "promises"
))
```

### Bioconductor packages (for pathway analysis)

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "clusterProfiler",
  "org.Mm.eg.db",   # mouse gene annotations
  "org.Hs.eg.db"    # human gene annotations
))
```

### Optional packages

These extend the app's functionality but are not required to run the core pipeline:

| Package | Feature | Install |
|---------|---------|---------|
| `BPCells` | On-disk storage for large datasets (>150K spots) | `remotes::install_github("bnprks/BPCells/r")` — requires `libhdf5-dev` |
| `Azimuth` | Cell type prediction against reference atlases | `remotes::install_github("satijalab/azimuth")` |
| `Banksy` | Spatially-aware clustering | `BiocManager::install("Banksy")` |
| `SeuratWrappers` | BANKSY–Seurat interface | `remotes::install_github("satijalab/seurat-wrappers")` |
| `harmony` | Harmony batch correction (multi-sample) | `install.packages("harmony")` |

---

## Installation

```bash
git clone https://github.com/<your-username>/<your-repo>.git
cd <your-repo>
```

No build step is required. The app is launched directly from R.

---

## Running the App

The app requires two processes running simultaneously: the **Shiny app** and the **background worker**. The worker handles all computationally expensive steps (normalisation, PCA, clustering, marker finding) so the UI stays responsive.

### Step 1 — Start the worker

Open a terminal and run:

```bash
Rscript worker.R /path/to/queue /path/to/results
```

The worker will print `Lurking for jobs` when ready. Keep this terminal open — the worker must be running before you submit any analysis jobs from the app.

The queue and results directories will be created automatically if they do not exist. On an HPC system you can run the worker as a background job:

```bash
nohup Rscript worker.R /path/to/queue /path/to/results > worker.log 2>&1 &
```

The paths can also be set via environment variables instead of command line arguments:

```bash
export VISIUMHD_QUEUE_DIR=/path/to/queue
export VISIUMHD_RESULTS_DIR=/path/to/results
Rscript worker.R
```

### Step 2 — Launch the Shiny app

In a second terminal (or in RStudio):

```r
shiny::runApp("/path/to/repo")
```

Or set the queue and results paths explicitly if they differ from the defaults:

```r
Sys.setenv(VISIUMHD_QUEUE_DIR   = "/path/to/queue")
Sys.setenv(VISIUMHD_RESULTS_DIR = "/path/to/results")
shiny::runApp("/path/to/repo")
```

The app will open in your browser at `http://127.0.0.1:<port>`.

---

## Using the App

### Option A — Load the example dataset
Click **Load Example Dataset** on the landing page to load a pre-processed mouse brain Visium HD object. No upload is required. This is the fastest way to explore the app's features.

### Option B — Upload your own data
Click **Upload Your Own File** and select either:
- A **zipped Spaceranger output folder** (containing `spatial/` and `raw_feature_bc_matrix/` or `filtered_feature_bc_matrix/`) — the app will build the Seurat object automatically
- A **pre-processed Seurat `.rds` file** — must be unclustered and unprocessed (no existing PCA/UMAP)

Select the bin size (8µm or 16µm) before uploading a Spaceranger zip. This setting has no effect when uploading an `.rds` file.

### Option C — Compare two samples
Click **Compare Two Samples** on the landing page to load two datasets for integrated multi-sample analysis.

### Analysis workflow

Once data is loaded, the app guides you through a two-step pipeline:

1. **QC & Setup tab** — review violin plots, adjust thresholds, select normalisation method, click **Filter & Run PCA**
2. **PCA & Clustering tab** — review the elbow plot, set number of PCs and clustering resolution, choose clustering method, click **Run UMAP & Clustering**
3. **Explore tab** — interact with spatial and UMAP plots; colour by cluster or gene expression; draw spatial selections
4. **Markers tab** — cluster markers table, heatmap, volcano plots, pathway analysis
5. **Cell Type Labelling tab** — Azimuth annotation, PanglaoDB lookup, custom cluster renaming
6. **Export** — download figures, processed object, marker tables, and pipeline settings

Progress is shown in the **Pipeline Log** at the bottom of the sidebar. A job status card appears at the top of the sidebar during background computation.

---

## HPC / Server Deployment

On an HPC cluster, it is recommended to:

1. Run the worker as a SLURM job, setting `--cpus-per-task` to the desired parallelism — the worker reads `SLURM_CPUS_PER_TASK` automatically
2. Use a shared filesystem path for the queue and results directories so both the Shiny server process and the worker process can access them
3. Run the Shiny app via `shiny-server` or similar, pointing at the repository directory

For multi-user deployments, each user session uses a unique `session_id` to namespace its files within the results directory, so concurrent users do not interfere with each other.

---

## Notes on Large Datasets

- Datasets with **>150K spots**: counts are automatically converted to BPCells on-disk format before normalisation if BPCells is installed. All spots are retained.
- Datasets with **>300K spots**: sketch-based clustering is automatically enabled (50K representative spots clustered, labels projected back to all spots).
- SCTransform is **not compatible** with BPCells on-disk storage and will be automatically switched to NormalizeData for large datasets. This is not a limitation — NormalizeData is the recommended method for Visium HD spatial data.

---

## Citation

If you use this app in your research, please cite the underlying tools:

- **Seurat v5**: Hao et al. (2024) *Nature Methods*
- **BANKSY**: Singhal et al. (2024) *Nature Genetics*
- **Azimuth**: Hao et al. (2021) *Cell*
- **BPCells**: Corces et al. (2024)
- **clusterProfiler**: Wu et al. (2021) *The Innovation*
- **Harmony**: Korsunsky et al. (2019) *Nature Methods*

---

## Running the pipeline

On a SLURM cluster:

```bash
sbatch submit.sh
```

Or directly with Nextflow:

```bash
nextflow run main.nf -profile slurm -resume
```

Input paths (data directory, output directory, container, library paths) are set
in `nextflow.config`; edit these to match your environment before running.

---

## Data availability

This project uses publicly available Visium HD colorectal cancer data from
Oliveira et al. (2025), *Nature Genetics*. The raw spatial data and the RCTD
deconvolution results (`DeconvolutionResults_<SAMPLE>CRC.csv.gz`) are not
redistributed here; obtain them from the original publication's data release.

---

## Environment & Reproducibility

Full `sessionInfo()` output captured at the end of each run is provided in `env/`,
alongside a fixed random seed (`set.seed(42)`) and the fixed Bioconductor 3.19
container, recording the exact computational environment behind the reported results.

---

## Citation

If you use this pipeline, please cite the underlying tools:

- **Seurat v5**: Hao et al. (2024) *Nature Methods*
- **BANKSY**: Singhal et al. (2024) *Nature Genetics*
- **RCTD (spacexr)**: Cable et al. (2022) *Nature Biotechnology*
- **clusterProfiler**: Wu et al. (2021) *The Innovation*
- **Source dataset**: Oliveira et al. (2025) *Nature Genetics*

---

## License

See `LICENSE`.

---

## Author

Michael Tuft — MSc Applied Bioinformatics, King's College London, in partnership
with the Genomics R&D Platform, Guy's Hospital, NHS
