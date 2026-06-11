# CosMx Spatial Transcriptomics
Analysis repository for **Nanostring CosMx SMI** spatial transcriptomics data

---

## Overview

This repository contains the full computational and statistical analysis pipeline for CosMx spatial transcriptomics data.

The work covers the complete analysis workflow from raw flat-file ingestion through quality control, dimensionality reduction, clustering, and cell-type annotation, with an emphasis on spatially-aware methods suited to multiplexed imaging transcriptomics.

---

## Repository Structure

```
.
├── main.Rmd                        # Primary analysis notebook
├── src/
│   ├── custom_functions.R          # Franz's Seurat-native CosMx pipeline
│   ├── prep_cosmx.R                # Giuseppe's reference data prep pipeline
│   └── gb_model.Rmd                # Giuseppe's full reference workflow notebook
├── data/
│   ├── flatFiles/                  # Raw CosMx flat files (gitignored)
│   │   └── TA649/
│   │       ├── TA649_exprMat_file.csv.gz
│   │       ├── TA649_metadata_file.csv.gz
│   │       ├── TA649_tx_file.csv.gz
│   │       ├── TA649-polygons.csv.gz
│   │       ├── TA649_fov_positions_file.csv.gz
│   │       └── TA-649 Autoimmune Hepatitis biopsy.xlsx  # TMA construction key
│   ├── seurat_obj.RDS              # QC-filtered Seurat snapshot (Git LFS, 335 MB)
│   └── md5sum/
│       └── md5sum_flatFiles.csv    # Integrity checksums for raw flat files
├── doc/
│   └── session_custom_functions.md # Notes on custom function development
├── .gitattributes                  # Git LFS tracking rules
├── .gitignore
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
- Known AIH/liver marker panels
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
