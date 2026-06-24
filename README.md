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
├── README.md

├── LICENSE
├── main.nf                              # Nextflow pipeline: runs both arms per sample

├── nextflow.config                      # Profiles (SLURM/local), container, resources

├── submit.sh                            # SLURM submission wrapper

├── visiumhd_reproducibility_pipeline.Rmd  # Reproducibility validation

├── scripts/

│   ├── phase1_cluster.R                 # Standard (expression-only) clustering arm

│   ├── banksy_cluster.R                 # Spatially-aware (BANKSY) clustering arm

│   ├── 00_helpers.R                     # Shared helper functions

│   ├── 01_cross_sample_overview.R       # Cross-sample concordance / coherence

│   ├── 02_highlight_cluster_analysis.R  # (verify name + purpose)

│   └── 03_target_cluster.R              # (verify name + purpose)

├── shiny_app/                           # Redacted Shiny app

│   ├── analysis_module.R                # Spatial/UMAP plots, marker analysis, volcano plots

│   ├── annotation.R                     # PanglaoDB cell type lookup

│   ├── app.R                            # Sources all scripts to run the app

│   ├── azimuth_module.R                 # Azimuth cell type prediction

│   ├── export_module.R                  # Figure export, methods paragraph, save object

│   ├── global.R                         # Shared constants, helper functions, CSS 

│   ├── help_tab.R                       # Provides information for users

│   ├── multi_sample_module.R            # Multi-sample loading and comparison 

│   ├── qc_module.R                      # QC filtering, normalisation, PCA

│   ├── queue_module.R                   # Job queue polling and status management

│   ├── ui_server.R                      # App UI layout and main server function

│   └── worker.R                         # Background compute worker (run independently)

└── env/                                 # sessionInfo() records (add to back §ref in report)

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
each run is in `env/`. Nextflow ≥ 23.10.0 is required.

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
