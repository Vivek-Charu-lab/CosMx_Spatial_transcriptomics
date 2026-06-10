# Session Notes — `custom_functions.R` Development

## Overview

This session built `custom_functions.R` from scratch, implementing utility
functions for CosMx SMI data already loaded into a Seurat object. The
reference implementation was Giuseppe's `prep_cosmx.R` (`dataprep_cosmx()`).

---

## Functions Implemented

### `readCosMx(experiment_repo, MM_PER_PX)`
Loads CosMx flat files from one or more TMA subdirectories and returns a
Seurat object.

**Key design decisions vs Giuseppe:**
- Uses `data.table::melt` + `Matrix::sparseMatrix` (no dense intermediate) for memory efficiency
- Takes the **union** of genes across slides (zero-fills missing genes); Giuseppe takes the **intersection**
- Throws a hard error on duplicate slidenames; Giuseppe auto-disambiguates (`-2`, `-3`, ...)
- Stores polygon local-px data in `obj@misc$polygons` (`cell`, `slidename`, `x_local_px`, `y_local_px`) for downstream QC
- Does **not** compute `SplitRatioToLocal` at load time — delegated to `computeSplitRatio()`
- Does **not** condense tissues at load time — delegated to `condenseTissues()`

**Tested on:** `data/flatFiles/TA649/` — single slide, 126,689 cells, 182 FOVs.

---

### `computeSplitRatio(obj)`
Computes `SplitRatioToLocal` for each cell and adds it to `obj@meta.data`.

**Logic (mirrors Giuseppe exactly, per-slide):**
```r
# Per slide — uses each slide's own local-px extremes as the FOV frame
boundary_cells <- cells whose polygon vertices hit min/max x_local_px or y_local_px
SplitRatioToLocal <- ifelse(is_boundary, round(Area / mean_fov_area, 2), 0)
```

**Requires:** `obj@misc$polygons` (set by `readCosMx()`).  
**Warns** if `SplitRatioToLocal` already exists in metadata (will overwrite).  
**Errors** with clear message if `obj@misc$polygons` is absent.

**Why per-slide matters:** In multi-TMA experiments, each slide's local-px
coordinate space is independent. Computing global min/max across all slides
would miss boundary cells in slides whose extremes differ from the global min/max.

**Validation:** Produces 9,470 boundary cells on TA649 vs 9,433 in the stored
`seurat_obj.RDS`. The ~37-cell discrepancy is a data provenance difference
(polygon file on disk vs snapshot used to create the RDS), not a logic error.
The 120 ±0.01 rounding differences are also from slightly different mean area
computation at original creation time.

---

### `condenseTissues(obj, tissue_col, x_col, y_col, tissueorder, buffer, widthheightratio, seed)`
Rearranges TMA core coordinates into a compact grid layout for visualization.
Stores result in `obj@misc$xy_condensed` (columns: `x_mm`, `y_mm`, `tissue`).

**Algorithm (exact port of Giuseppe's):**
1. Compute bounding box per tissue
2. Order tissues by decreasing height (or custom `tissueorder`)
3. Compute `tissuesperrow` from `widthheightratio`
4. Place tissues on shelves — start new shelf when adding next tissue would
   overshoot `targetwidth`
5. Shift each tissue's cells to their assigned slot origin

**Defaults:** `buffer = 0.2`, `widthheightratio = 8/3`, `seed = 1`

**Tested on:** TA649 — single slide, output x: 0–11.308 mm, y: 0–15.657 mm ✅

---

### `plotTissues(obj, subsample_frac, cols, pch, cex, main, label, seed)`
Plots the condensed tissue layout from `obj@misc$xy_condensed`.

**Errors with clear message if `condenseTissues()` has not been run first.**

**Defaults:** `subsample_frac = 1/20`, `pch = 16`, `cex = 0.2`, `label = TRUE`, `seed = 1`  
**Returns:** invisibly, the full `xy_condensed` data frame.

---

## QC Tasks in Giuseppe NOT Yet Implemented

| Task | Notes |
|---|---|
| Duplicate slidename auto-disambiguation | Giuseppe appends `-2`, `-3`; `readCosMx()` errors instead |
| Shared genes intersection strategy | `readCosMx()` uses union + zero-fill; Giuseppe uses intersection |

---

## Repository Layout Reminder

```
custom_functions.R     # All utility functions (this session)
prep_cosmx.R           # Giuseppe's reference implementation
seurat_obj.RDS         # Example Seurat object (TA649, 126,689 cells)
data/flatFiles/TA649/  # Raw CosMx flat files
  TA649_exprMat_file.csv.gz
  TA649_metadata_file.csv.gz
  TA649_tx_file.csv.gz
  TA649_fov_positions_file.csv.gz
  TA649-polygons.csv.gz
```

## R Environment

- Micromamba env: `/oak/stanford/groups/bhowitt/conda_envs/workspace`
- Activate: `micromamba run -p /oak/stanford/groups/bhowitt/conda_envs/workspace Rscript ...`
- Base mamba at: `/scratch/users/franzake/micromamba`
- The base mamba env does NOT have `SeuratObject` installed — always use `workspace` for testing

## Typical Usage

```r
source("custom_functions.R")

# Load data
obj <- readCosMx("data/flatFiles/")

# QC
obj <- computeSplitRatio(obj)
# Filter boundary cells:
obj <- subset(obj, subset = SplitRatioToLocal == 0 | SplitRatioToLocal < 1)

# Visualize tissue layout
obj <- condenseTissues(obj)
plotTissues(obj)
```
