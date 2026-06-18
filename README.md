# Visium HD Spatial Clustering: Comparative Pipeline & Reproducibility

This repository contains the reproducible analysis pipeline for an MSc research
project comparing standard (expression-only) and spatially-aware (BANKSY)
clustering of 10x Genomics Visium HD spatial transcriptomics data, applied to
three colorectal cancer (CRC) samples.

The project has two components:

1. **A comparative clustering pipeline** (this repository) — a Nextflow workflow
   that processes each sample through two matched clustering arms and a set of
   downstream comparison analyses, plus a reproducibility-validation pipeline.
2. **An interactive R Shiny application** for end-to-end Visium HD analysis. The
   application's full source is retained as proprietary and is **not** published
   here; redacted excerpts are provided separately for assessment purposes only
   (see *Scope* below).

---

## Scope: what is and isn't in this repository

**Included (fully reproducible):**
- The Nextflow comparative-clustering pipeline (`main.nf`, `nextflow.config`, `submit.sh`)
- The per-sample clustering scripts for both arms (`scripts/phase1_cluster.R`, `scripts/banksy_cluster.R`)
- The cross-sample comparison and deep-dive analysis scripts (`scripts/00`–`03`)
- The reproducibility-validation pipeline (`visiumhd_reproducibility_pipeline.Rmd`)

**Not included:**
- The interactive Shiny application's source code, which is proprietary. Redacted
  excerpts demonstrating the implementation are provided to examiners for
  assessment only and are not runnable.
- The raw Visium HD data and deconvolution results — these are third-party and
  available from the original publication (see *Data Availability*).

---

## Repository structure
├── README.md

├── LICENSE

├── main.nf                              # Nextflow pipeline: runs both arms per sample

├── nextflow.config                      # Profiles (SLURM/local), container, resources

├── submit.sh                            # SLURM submission wrapper

├── scripts/

│   ├── phase1_cluster.R                 # Standard (expression-only) clustering arm

│   ├── banksy_cluster.R                 # Spatially-aware (BANKSY) clustering arm

│   ├── 00_helpers.R                     # Shared helper functions

│   ├── 01_cross_sample_overview.R       # Cross-sample concordance / coherence overview

│   ├── 02_highsplit_cluster_analysis.R  # High-complexity cluster analysis

│   └── 03_target_cluster_deep_dive.R    # Per-sample flagship-cluster deep dive

└── visiumhd_reproducibility_pipeline.Rmd  # Reproducibility validation (clustering, RCTD, Azimuth, BANKSY rotation)

*(Full `sessionInfo()` records for each run are provided in `env/` — see
**Environment & Reproducibility**.)*

---

## The analysis

Each sample (P1, P2, P5) is processed independently through two clustering arms
that share identical preprocessing, so that any difference in the resulting
clusters is attributable to the spatial augmentation rather than to differing
settings:

- **Standard arm** (`phase1_cluster.R`): expression-only clustering. Spots are
  filtered (≥100 UMIs, ≥50 features), log-normalised (scale factor 10,000), 3,000
  variable features, 20 principal components, shared-nearest-neighbour graph,
  Louvain clustering at resolution 0.8.
- **BANKSY arm** (`banksy_cluster.R`): identical preprocessing, with the
  expression matrix augmented by neighbourhood features (BANKSY; λ = 0.2,
  k_geom = 50, `kNN_median`) before PCA and the same SNN + Louvain clustering.

Cell-type composition is annotated from the published RCTD deconvolution results
released with the source dataset (matched per barcode). The comparison scripts
(`00`–`03`) quantify cross-arm concordance, spatial coherence, and per-cluster
composition.

The QC thresholds (≥100 UMIs, ≥50 features) are matched to the source paper's
deconvolution threshold for comparability.

---

## Requirements

Analyses were run in **R 4.4.1** inside the standard **Bioconductor 3.19
Singularity image** (`bioc_3.19.sif`, Ubuntu 22.04.4 LTS). Key packages: Seurat
v5, BANKSY, SeuratWrappers, spacexr (RCTD), Azimuth, clusterProfiler. The
complete `sessionInfo()` for each run is in `env/`.

Nextflow ≥ 23.10.0 is required to run the pipeline.

---

## Running the pipeline

The pipeline processes all three samples through both arms in parallel.

On a SLURM cluster:

```bash
sbatch submit.sh
```

Or directly with Nextflow:

```bash
nextflow run main.nf -profile slurm -resume
```

For local execution (single machine, no scheduler):

```bash
nextflow run main.nf -profile standard
```

Input paths (data directory, output directory, container, and library paths) are
set in `nextflow.config`. Edit these to match your environment before running.

---

## Data availability

This project uses publicly available Visium HD colorectal cancer data from:

> Oliveira et al. (2025), *Nature Genetics*.

The raw spatial data and the RCTD deconvolution results
(`DeconvolutionResults_<SAMPLE>CRC.csv.gz`) are not redistributed here and should
be obtained from the original publication's data release. Place them in the
locations expected by `nextflow.config` before running.

---

## Environment & Reproducibility

Full `sessionInfo()` output captured at the end of each run is provided in `env/`:

env/

├── session_info_phase1.txt          # Standard-arm runs

├── session_info_banksy.txt          # BANKSY-arm runs

└── session_info_reproducibility.txt # Reproducibility-validation pipeline
Together with the fixed container (`bioc_3.19.sif`) and a fixed random seed
(`set.seed(42)`), these record the exact computational environment behind the
reported results.

---

## Citation

If you use this pipeline, please cite the underlying tools:

- **Seurat v5**: Hao et al. (2024) *Nature Methods*
- **BANKSY**: Singhal et al. (2024) *Nature Genetics*
- **RCTD (spacexr)**: Cable et al. (2022) *Nature Biotechnology*
- **Azimuth**: Hao et al. (2021) *Cell*
- **clusterProfiler**: Wu et al. (2021) *The Innovation*
- **Source dataset**: Oliveira et al. (2025) *Nature Genetics*

---

## License

See `LICENSE`.

---

## Author

Michael Tuft — MSc Applied Bioinformatics, King's College London, in partnership
with the Genomics Research Platform, Guy's Hospital.
