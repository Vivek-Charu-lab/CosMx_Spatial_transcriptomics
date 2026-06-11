# Copilot Instructions

## Project Overview

Spatial transcriptomics analysis lab (Vivek Lab). Contains raw **Nanostring CosMx SMI** data for sample **TA649** (Human RNA 6k Discovery panel, ~6,000 genes; 126,689 cells, 182 FOVs). Primary analysis is in R, run interactively in Posit/RStudio.

## Running the Analysis

No formal build/test system. Analysis is executed by knitting R Markdown notebooks:

```r
# From R console (project root):
rmarkdown::render("main.Rmd")      # Main analysis notebook (uses custom_functions.R)
rmarkdown::render("gb_model.Rmd")  # Giuseppe's reference pipeline

# Or use RStudio/Posit's Knit button
```

Parallelism is configured at the top of `main.Rmd`:
```r
plan("multisession", workers = 4)
options(future.globals.maxSize = 20 * 1024^3)  # 20 GB
```

## R Environment

- **Required conda env**: `/oak/stanford/groups/bhowitt/conda_envs/workspace` (micromamba)
- Activate for CLI use: `micromamba run -p /oak/stanford/groups/bhowitt/conda_envs/workspace Rscript ...`
- The base mamba env does **not** have `SeuratObject` installed — always use `workspace` for R work

## Architecture

### Two parallel pipelines

**1. `custom_functions.R`** — Franz's Seurat-native pipeline (primary for `main.Rmd`):

- `readCosMx(experiment_repo)` → **Seurat object** with:
  - Assay `RNA`: gene counts (genes × cells); `negprobes` and `falsecode` as separate assays
  - One spatial FOV per slide named by `slidename` (e.g., `obj[["TA649"]]`), with centroids and segmentation polygons in mm
  - `obj@misc$polygons`: polygon local-px coords for downstream QC (`cell`, `slidename`, `x_local_px`, `y_local_px`)
- `computeSplitRatio(obj)` — adds `SplitRatioToLocal` to `@meta.data`; requires `obj@misc$polygons`; cells with value `> 1` are strong filter candidates
- `condenseTissues(obj)` — packs TMA cores into a compact grid; stores result in `obj@misc$xy_condensed` (`x_mm`, `y_mm`, `tissue`); **called automatically inside `readCosMx()`**
- `TmaPlot(obj, col.by, pt.size, main, dark, cols, facet, subsample_frac, legend.max.levels, seed)` — primary spatial viz; plots condensed mm coords coloured by any metadata column, faceted by tissue; auto-calls `condenseTissues()` if `obj@misc$xy_condensed` is absent; continuous columns use viridis, discrete use categorical palette; legend hidden when `> legend.max.levels` (default 50) unique values
- `computeSBR(obj, bandwidth, weight_cutoff)` — adds `SBR` and `log2SBR` to `@meta.data`; Gaussian-kernel spatial smoothing of `nCount_RNA` / `nCount_negprobes` ratio; uses **condensed coordinates** (`obj@misc$xy_condensed`); requires `dbscan` package
- `computeFOVintegrity(obj, panel_name, max_prop_loss, max_totalcounts_loss)` — wraps `gbspatial:::runFOVQC()`; adds `flag_fov_integrity` (logical) and `fov_signal_loss` (log2 FC) to `@meta.data`; uses **condensed coordinates** (`obj@misc$xy_condensed`); stores full result in `obj@misc$fov_integrity`; default panel `"Hs_6k"`
- `plotFOVintegrity(obj)` — renders two base-graphics plots from `obj@misc$fov_integrity` via internal `gbspatial` functions: flagged FOV map and signal-loss spatial plot; returns invisible list of recorded plots
- `QCdistributions(obj, count_min/max, area_min/max, sbr_min/max, split_ratio_min/max, bins, title)` — histogram panel for nCount_RNA, Area, log2SBR, SplitRatioToLocal; returns named list of ggplots
- `QCplots(obj, count_min/max, area_min/max, sbr_min/max, split_ratio_min/max, title)` — scatter-plot panel (nCount vs nFeature, nCount vs Area, SBR vs nCount, SplitRatio vs Area) with optional red threshold lines; returns named list of ggplots
- `QCfiltering(obj, count_min, count_max, area_min, area_max, split_ratio_min, split_ratio_max, sbr_min, sbr_max, panel_name, fov_integrity_threshold, do_*)` — full QC pipeline; adds per-cell flag columns (`flag_nCount`, `flag_area`, `flag_fov_integrity`, `flag_boundary`, `flag_sbr`, `flag_overall`) and optionally filters; individual steps can be toggled with `do_*` arguments
- `parseTMAmap(path, sheet, first_col, last_col, rows_per_core, first_pid_row, flip_rows, flip_cols)` — reads a TMA construction-key Excel file (default sheet `"TA649 constructionKey"`) and returns a character matrix of patient IDs (rows = TMA rows, cols = TMA cols); strips everything after the first comma in each cell (handles multi-value entries); supports flipping rows/cols for orientation correction
- `assignCoresToCells(obj, tma_map, fov_positions, fov_size, core_drift_tolerance, min_cells, sample_id_col)` — wraps `gbspatial::assign_fovs_to_cores()`; maps each FOV to a TMA core grid position and propagates to cells; adds `core_str` (e.g. `"C3R2"`), `core_col`, `core_row`, and `sample_id_col` (default `"sample_id"`) to `@meta.data`; unassigned cells receive `NA`; full FOV-level result stored in `obj@misc$fov_core_mapping`; diagnostic plots in `obj@misc$fov_core_plots`

**2. `prep_cosmx.R`** — Giuseppe's reference pipeline (see also `gb_model.Rmd`):

- `dataprep_cosmx(myflatfiledir)` → **list** with `counts`, `negcounts`, `falsecounts`, `metadata`, `xy`, `polygons`
- Reads counts in dense chunks (memory-safe but slow); `readCosMx()` uses melt + `sparseMatrix` directly — no dense intermediate
- `gb_model.Rmd` shows the full downstream workflow: QC via `gbspatial::run_spatial_qc()` → UMAP+clustering via `scPearsonPCA` → cell typing via `HieraType` → spatial polygon plots

### Key design differences: `readCosMx` vs `dataprep_cosmx`

| Behavior | `readCosMx` | `dataprep_cosmx` |
|---|---|---|
| Output type | Seurat object | List |
| Multi-slide gene handling | Union + zero-fill | Intersection |
| Duplicate slidenames | Hard error | Auto-disambiguate (`-2`, `-3`) |
| `SplitRatioToLocal` | Computed separately via `computeSplitRatio()` | Computed inline at load time |
| Tissue condensation | Stored in `obj@misc$xy_condensed` | Returns condensed `xy` matrix |

### Data flow

```
data/flatFiles/TA649/   →   readCosMx()           →   Seurat obj   →   main.Rmd
                        →   dataprep_cosmx()       →   list obj     →   gb_model.Rmd
data/seuratObject_TA.649.RDS  ← pre-built (2.6 GB), skip ingestion for quick work
seurat_obj.RDS          ← another pre-built snapshot at repo root
```

## Typical Usage

```r
source("custom_functions.R")

obj <- readCosMx("data/flatFiles/")
# condenseTissues() runs automatically inside readCosMx()

# Visualize TMA layout
TmaPlot(obj)                             # colour by orig.ident (default)
TmaPlot(obj, col.by = "nCount_RNA")      # colour by continuous metadata

# Full QC pipeline (mirrors Giuseppe's run_spatial_qc())
obj <- computeSplitRatio(obj)
obj <- computeSBR(obj)
obj <- computeFOVintegrity(obj)
plotFOVintegrity(obj)                    # inspect flagged FOVs spatially

# Inspect with plots
QCdistributions(obj)
QCplots(obj, count_min = 30, area_max = 70e3, sbr_min = 0, sbr_max = 10)

# Filter (sets flag_overall; removes flagged cells when filter = TRUE)
obj_qc <- QCfiltering(obj,
  count_min              = 30,
  area_max               = 70e3,
  sbr_min                = 0,
  sbr_max                = 10,
  fov_integrity_threshold = 0.6
)

# Assign TMA core identities — parse TMA map from Excel first
tma_map <- parseTMAmap("data/flatFiles/TA649/TA-649 Autoimmune Hepatitis biopsy.xlsx")
fov_pos <- "data/flatFiles/TA649/TA649_fov_positions_file.csv.gz"
obj_qc <- assignCoresToCells(obj_qc, tma_map = tma_map, fov_positions = fov_pos)
# Adds core_str, core_col, core_row, sample_id to @meta.data
TmaPlot(obj_qc, col.by = "sample_id")
```

## Downstream Analysis Pattern (`main.Rmd`)

After QC and core assignment, `main.Rmd` follows a standard Seurat workflow:

```r
# Parallelism (configured once at top of notebook)
plan("multicore", workers = 6)
options(future.globals.maxSize = 60 * 1024^3)

# Dimensionality reduction
obj_qc <- FindVariableFeatures(obj_qc, nfeatures = 2000)
obj_qc <- SCTransform(obj_qc, assay = "RNA",
                      vars.to.regress = c("nCount_RNA", "nFeature_RNA"),
                      clip.range = c(-10, 10))
obj_qc <- RunPCA(obj_qc, npcs = 30, assay = "SCT")
obj_qc <- RunUMAP(obj_qc, dims = 1:NPCS)   # NPCS chosen from ElbowPlot

# Clustering (sweep resolutions, inspect with clustree)
obj_qc <- FindNeighbors(obj_qc, dims = 1:NPCS)
obj_qc <- FindClusters(obj_qc, resolution = seq(0, 1, 0.1))
clustree(obj_qc)
obj_qc$seurat_clusters <- obj_qc$SCT_snn_res.0.5   # set chosen resolution

# Cell-type annotation
markers <- data.table(FindAllMarkers(obj_qc, assay = "SCT", only.pos = TRUE))
```

Known cell-type marker panels for this dataset (Autoimmune Hepatitis / liver):
- Kupffer/Macrophages: `C1QC`, `C1QB`, `C1QA`, `MARCO`, `CD163`
- Hepatocytes: `SLC27A5`, `ADH4`, `SLC2A2`, `ARG1`, `CYP4A11`
- HSCs: `COL1A1`, `COL3A1`, `DCN`, `TAGLN`, `SPARC`
- Cholangiocytes: `KRT7`, `KRT8`, `TACSTD2`, `MMP7`, `SPP1`
- Pericentral hepatocytes: `GLUL`, `TDO2`, `FKBP5`; Periportal: `SAA1`, `CRP`, `FBP1`
- Plasma cells: `IGKC`, `IGHG1/2`, `JCHAIN`, `MZB1`, `IGHA1`

## Key Conventions

### Cell & FOV identifiers
- **Global cell ID**: `paste0("c_", slide_index, "_", fov, "_", cell_ID)` — `slide_index` is the 1-based position from directory enumeration
- **Global FOV ID**: `paste0("s", slide_index, "f", fov)` — e.g., `s1f42`
- `cell_ID` (integer) is only unique within a FOV — always use the string global cell ID for cross-file joins

### Coordinate systems
- `*_local_px`: within-FOV pixels
- `*_global_px`: slide-level pixels
- `*_slide_mm` / `x_mm` / `y_mm`: millimeters (conversion: `px × 0.120280945 / 1000`)
- Seurat spatial FOVs use mm coordinates; `obj@misc$polygons` keeps local-px for `computeSplitRatio()`

### Counts orientation
- Flat files: cells × genes (rows = cells)
- Seurat expects genes × cells — always `t(counts)` when calling `CreateSeuratObject()`

### Control probe separation
- `grepl("Negative", colnames(counts))` → `negprobes` assay
- `grepl("SystemControl", colnames(counts))` → `falsecode` assay
- Never include these in gene expression analyses

### `SplitRatioToLocal` boundary detection
- Computed **per slide** using each slide's own `x_local_px`/`y_local_px` extremes — critical in multi-TMA experiments where each slide's coordinate space is independent
- Validation on TA649: 9,470 boundary cells detected vs 9,433 in the pre-built RDS — ~37-cell discrepancy is a data provenance difference (polygon file on disk vs the snapshot used to build the RDS), not a logic error

### `TmaPlot` axis convention
- X and Y are **intentionally swapped** in the ggplot call: `aes(x = y_mm, y = x_mm)` — this matches the physical orientation of CosMx TMA slides

### Data loading
- All flat files are gzip-compressed CSVs; `data.table::fread()` handles `.gz` automatically
- `fov_positions_file.csv.gz` contains two concatenated tables — column headers reset mid-file, parse carefully
- Verify integrity with `data/md5sum/md5sum_flatFiles.csv` before analysis

## CosMx Flat-File Schema

### `*_exprMat_file.csv.gz`
- Columns: `fov`, `cell_ID`, then one column per gene (integer UMI counts)

### `*_metadata_file.csv.gz`
- Cell identifiers: `fov`, `cell_ID` (integer), `cell_id` / `cell` (string `c_<slide>_<fov>_<cellID>`)
- Morphology: `Area`, `AspectRatio`, `Circularity`, `Eccentricity`, `Perimeter`, `Solidity`, `NucArea`, `NucAspectRatio`
- Fluorescence: `Mean.DAPI`/`Max.DAPI` (nuclei), `Mean.PanCK`/`Max.PanCK` (epithelial), `Mean.CD45`/`Max.CD45` (immune), `Mean.Membrane`/`Max.Membrane`, `Mean.G`/`Max.G`
- RNA QC: `nCount_RNA`, `nFeature_RNA`, `unassignedTranscripts`, `median_RNA`, `RNA_quantile_*`
- Negative probe QC: `nCount_negprobes`, `nFeature_negprobes`, `median_negprobes`
- False code QC: `nCount_falsecode`, `nFeature_falsecode`, `median_falsecode`
- Spatial: `CenterX_global_px`, `CenterY_global_px`, `CenterX_local_px`, `CenterY_local_px`, `Area.um2`

### `*_tx_file.csv.gz`
- Per-transcript: `fov`, `cell_ID`, `cell`, `x_local_px`, `y_local_px`, `x_global_px`, `y_global_px`, `z`, `target` (gene), `CellComp` (subcellular compartment)
- `cell_ID = 0` → unassigned transcript

### `*-polygons.csv.gz`
- Polygon vertices: `fov`, `cellID`, `cell`, `x_local_px`, `y_local_px`, `x_global_px`, `y_global_px`

## Ecosystem Context

- **Core R packages**: `Seurat` (v5), `SeuratObject`, `Matrix`, `data.table`, `future`/`future.apply`, `ggplot2`, `dplyr`, `dbscan` (required by `computeSBR()`), `openxlsx` (required by `parseTMAmap()`), `clustree` (resolution sweep visualization)
- **Lab-specific packages**: `gbspatial` (Giuseppe's — `dataprep_cosmx()`, `run_spatial_qc()`, `assign_fovs_to_cores()`), `scPearsonPCA` (quasi-Poisson PCA + UMAP), `HieraType` (cluster markers + cell typing)
- **Python alternative**: `squidpy`, `scanpy`, `spatialdata` with `spatialdata-io` CosMx reader

## Session Log

### 2026-06-10 / 2026-06-11

**Completed this session:**
- `main.Rmd` restructured: 6 numbered `##` sections, `###` subsections, consistent `<-` assignment, horizontal rules
- Added **Ideas & Tasks to Develop** section (top of notebook) with 3 items from team meeting:
  1. Transfer learning from Nanostring liver reference (Azimuth / `FindTransferAnchors`)
  2. Manual annotation + CellMarker 2.0 cross-reference
  3. `plotSegmentation()` function using polygon vertices + ggplot2 — discuss design with Karan first
- Bug fix: added missing `computeSplitRatio(usc)` call in §4.1
- Removed orphaned `gbspatial::dataprep_cosmx` scratch chunk (broke knitting via `View()`)
- `README.md` created with project overview, pipeline stages, dependency tables, rclone/Box data access instructions
- `.gitignore` extended: `*.bak`, `*.orig`, notebook cache/HTML/PDF, `*.rds`
- Repo transferred: `franzake/` → `Vivek-Charu-lab/CosMx_Spatial_transcriptomics`
- Git remote updated locally; `~/.gitconfig` fixed (broken `gh` path: `miniforge3` → `/share/software/user/open/gh/2.88.1/bin/gh`)
- Git LFS configured: `data/seurat_obj.RDS` (335 MB) tracked and pushed
- Raw flat files (~2.6 GB) synced to Box: `CosMX AIH/TA649/flatFiles/`

**Current state of `main.Rmd`:**

| Section | Status |
|---|---|
| §1 Libraries & Configuration | ✅ |
| §2 Input Files | ✅ |
| §3 Data Loading | ✅ |
| §4.1 QC Metrics | ✅ |
| §4.2 FOV-to-Core Assignment | ✅ |
| §5.1–5.4 Dimensionality Reduction | ✅ (4 approaches compared) |
| §6.1 Marker Discovery | ✅ (needs cluster resolution choice) |
| §6.2 Cell-Type Annotation | 🔄 Feature plots for first 2 panels only |

**Next steps:**
- Choose final DR approach (§5.1–5.4) and run `FindClusters` + `clustree` resolution sweep
- Complete cell-type annotation for all 11 `jci_insight_markers` panels
- Transfer learning pipeline from Nanostring liver reference
- `plotSegmentation()` — discuss with Karan before implementing
- Manual annotation refinement with CellMarker 2.0

**Environment notes:**
- Git LFS binary must be in PATH for pushes: `export PATH="/share/software/user/open/git-lfs/2.4.0/bin:$PATH"`
- `gh` CLI: `/share/software/user/open/gh/2.88.1/bin/gh`
- HPC: Stanford Sherlock (`sh02-ln04`), project at `/scratch/users/franzake/CosMx_Spatial_transcriptomics`
