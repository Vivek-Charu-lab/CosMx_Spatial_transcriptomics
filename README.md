# CosMx Spatial Transcriptomics — TA649

Analysis repository for **Nanostring CosMx SMI** spatial transcriptomics data from the **Vivek Charu Lab** (Stanford University).

---

## Overview

This repository contains the full computational and statistical analysis pipeline for CosMx spatial transcriptomics data, currently applied to sample **TA649** — a TMA (Tissue Microarray) cohort of **Autoimmune Hepatitis (AIH)** liver biopsies profiled with the **Human RNA 6k Discovery panel** (~6,000 genes, 182 FOVs, ~126,000 cells).

The work covers the complete analysis workflow from raw flat-file ingestion through quality control, dimensionality reduction, clustering, and cell-type annotation, with an emphasis on spatially-aware methods suited to multiplexed imaging transcriptomics.

---

## Repository Structure

```
.
├── main.Rmd                  # Primary analysis notebook
├── src/
│   ├── custom_functions.R    # Franz's Seurat-native CosMx pipeline
│   └── prep_cosmx.R          # Giuseppe's reference pipeline (gbspatial)
├── data/
│   ├── flatFiles/            # Raw CosMx flat files (gitignored)
│   └── md5sum/               # Integrity checksums for raw data
├── doc/                      # Notes and documentation
└── CosMx_Spatial_transcriptomics.Rproj
```

---

## Analysis Pipeline

### 1. Data Ingestion
Raw CosMx flat files (expression matrix, metadata, polygon segmentations, FOV positions) are loaded into a **Seurat v5** object via `readCosMx()`. TMA core condensation is applied automatically.

### 2. Quality Control
Per-cell QC metrics computed and filtered:
- **nCount / nFeature RNA** — transcriptomic complexity
- **Cell Area** — segmentation size consistency
- **Signal-to-Background Ratio (SBR)** — Gaussian-kernel spatial smoothing of RNA vs negative probe counts
- **Split Ratio** — boundary fragmentation / over-segmentation detection
- **FOV Integrity** — field-of-view signal loss detection

### 3. TMA Core Assignment
FOVs are mapped to TMA core grid positions using the construction-key Excel file, propagating patient sample IDs (`sample_id`) to individual cells.

### 4. Dimensionality Reduction
Four embedding strategies are compared:
| Approach | Normalisation | Batch correction |
|---|---|---|
| Baseline | SCTransform | None |
| SCT + Harmony | SCTransform | Harmony |
| scPearsonPCA (native batch) | Quasi-Poisson PCA | Built-in |
| scPearsonPCA + Harmony | Quasi-Poisson PCA | Harmony |

### 5. Clustering & Cell-Type Annotation
Seurat graph-based clustering with resolution sweep (`clustree`). Cell types annotated using:
- Unsupervised marker discovery (`FindAllMarkers`)
- Known AIH/liver marker panels (Kupffer cells, hepatocytes, HSCs, cholangiocytes, plasma cells, etc.)
- CellMarker 2.0 database cross-reference
- Nanostring liver reference label transfer (in development)

---

## Dependencies

### R Packages
| Package | Role |
|---|---|
| `Seurat` v5 | Core single-cell / spatial framework |
| `scPearsonPCA` | Quasi-Poisson PCA for spatial count data |
| `gbspatial` | CosMx data prep & spatial QC (Vivek Lab) |
| `harmony` | Batch integration |
| `data.table`, `Matrix` | Fast data I/O |
| `ggplot2`, `patchwork` | Visualisation |
| `clustree` | Clustering resolution sweep |
| `dbscan` | Spatial smoothing (SBR computation) |
| `openxlsx` | TMA construction-key parsing |

### Environment
Analysis runs in a managed conda environment on the Stanford HPC (Sherlock):
```
/oak/stanford/groups/bhowitt/conda_envs/workspace
```

---

## Running the Analysis

Open the project in **RStudio / Posit** and knit `main.Rmd`, or from the R console:

```r
rmarkdown::render("main.Rmd")
```

Parallelism is configured at the top of `main.Rmd`:
```r
plan("multicore", workers = 6)
options(future.globals.maxSize = 60 * 1024^3)
```

---

## Author

**Franz Ake** — Vivek Charu Lab, Stanford University
