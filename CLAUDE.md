# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A **data-only workspace** on Sherlock (`$SCRATCH/COSMX_TA649`, also reachable as
`~/SCRATCH_REPO/COSMX_TA649` — same directory via symlink) holding the raw **NanoString CosMx SMI**
export for run `TA649`: an autoimmune-hepatitis (AIH) liver-biopsy TMA, Human RNA 6K Discovery panel
(`(1.1) Human RNA 6k Discovery; (1.0) 6K_Stanford_SB`), single slide `20260110_014400_S1`, 182 FOVs,
126,689 segmented cells.

There is no build system, package manifest, or git history here — `src/` holds plain scripts that
are sourced or submitted to Slurm, never run as a CLI.

## `src/` — the AnnData reader

- **`custom_funct.R`** — `readCosMxAnnData(experiment_repo, ...)`. Reads every immediate
  subdirectory of a CosMx repo as one TMA/slide and returns a Seurat v5 object. The expression
  matrix is streamed through Python `anndata` via `reticulate`; metadata and polygons are read in R
  with `data.table`. Verified **bit-identical** to the sibling `readCosMx()` (all three assays,
  nnz 125,813,845 / 272,849 / 821,901), so the sibling QC toolchain works on the result unchanged.
- **`test_custom_funct.R`** — diffs the two readers on the real slide (Part 1) and exercises the
  multi-TMA / zero-fill / `h5ad` / cache paths on a synthetic two-slide fixture (Part 2).
  `COSMX_TEST_PARTS=2` runs only the cheap part.
- **`run_test_custom_funct.sbatch`** — Slurm wrapper. Measured 16.1 GB peak, 5:14 wall.

Key facts for working on this code:

- **The R env ships no Python at all**, so reticulate silently falls back to `/usr/bin/python3`
  (no anndata) unless told otherwise. Pin
  `/scratch/users/franzake/micromamba/envs/despotx/bin/python` via `$COSMX_PYTHON`,
  `options(cosmx.python=)`, or the `python=` argument. `resolveCosMxPython()` fails loudly rather
  than letting this pass silently.
- The count matrix crosses the bridge as a **scipy CSC → `dgCMatrix`** with no dense step:
  `adata.X` is cells×genes CSR and `csr.T` is a zero-copy CSC, which is already Seurat's
  orientation, so nothing calls `Matrix::t()`. reticulate only registers converters for scipy
  sparse *matrix* classes — **not** the sparse *array* classes anndata ≥0.12 can return, which is
  why the despotx env (anndata 0.11.4) is pinned.
- `cache_dir=` drops a raw-counts `.h5ad` per slide, directly consumable by DeSpotX / ResolVI /
  scanpy. This is the main reason to prefer this reader over the pure-R one, which is ~2× faster
  (70 s vs 142 s) for an otherwise identical object.
- **`NegativeAdd` classification quirk:** the inherited `"Negative"` regex routes the feature named
  `NegativeAdd` into `negprobes`, but `plex-*.txt` classes it Endogenous and NanoString's own
  `nCount_RNA` counts it as RNA. Default split is therefore 6215/21/323, and negprobe QC metrics
  are slightly inflated. Kept as-is for parity with prior analyses; pass
  `control_patterns = c(negprobes = "^Negative[0-9]+$", falsecode = "SystemControl")` for the
  vendor's classification.
- Seurat **recomputes** `nCount_<assay>` whenever an assay is added, so comparing
  `obj$nCount_negprobes` against `colSums()` of that assay is tautological. Validate against the
  metadata CSV's pristine columns instead.

## Where the existing pipeline lives

Two sibling workspaces already analyze **this same TA649 dataset** and each carry their own
`CLAUDE.md` with far more detail. Read them before writing new analysis code here:

- `/scratch/users/franzake/CosMx_Spatial_transcriptomics/` — the R/Seurat-v5 pipeline
  (`main.Rmd` + `src/custom_functions.R`: `readCosMx`, QC metrics, `parseTMAmap`/`assignCoresToCells`,
  Harmony, clustering, annotation) plus the Python **DeSpotX** decontamination submodule and the
  RStudio-Server-on-Slurm launcher (`bin/launch_rstudio_server.sbatch`).
- `/scratch/users/franzake/Liver_CosMX/` — a later iteration (`main2.Rmd`) plus the **Minerva**
  image-build scripts that stitch the morphology TIFFs into an OME-TIFF and rasterize a cell-type
  label mask. Its `CLAUDE.md` documents the load-bearing coordinate/orientation conventions.

Do not duplicate `custom_functions.R` here — source it from the sibling workspace or copy
deliberately, and say which.

## Data layout

### `data/flatFiles/TA649/` — AtoMx flat-file export (2.6 GB, all `.csv.gz`)

| File | Contents |
|---|---|
| `TA649_exprMat_file.csv.gz` | cell × feature counts, 6561 cols = `fov,cell_ID` + 6559 features (142 MB) |
| `TA649_metadata_file.csv.gz` | per-cell morphology, IF intensities, QC quantiles, `CenterX/Y_global_px`, `cell_id` |
| `TA649_tx_file.csv.gz` | per-transcript records: `fov,cell_ID,cell,x/y_local_px,x/y_global_px,z,target,CellComp` (**2.6 GB gzipped — never load whole; stream it**) |
| `TA649-polygons.csv.gz` | segmentation polygon vertices, local **and** global px |
| `TA649_fov_positions_file.csv.gz` | FOV top-left global offsets, px and mm |
| `TA-649 Autoimmune Hepatitis biopsy.xlsx` | TMA key — sheets `List duplicated`, `TA649 constructionKey`, `TA649 PrintingKey (2)`; the input to `parseTMAmap()` (donor/tray/core → clinical annotation) |

Feature breakdown (from `plex-6atoipi5h2.txt`): **6216 Endogenous + 20 Negative + 323 SystemControl**.
`CellComp` levels in the tx file: `Nuclear`, `Membrane`, `Cytoplasm`, and empty (extracellular).
Cell IDs follow `c_<slide>_<fov>_<cell_ID>`, e.g. `c_1_1_1`.

### `data/DecodedFiles/TA649/20260110_014400_S1/` — instrument-level export (36 GB)

- `CellStatsDir/Morphology2D/` — per-FOV 5-channel morphology TIFFs,
  `20260110_014400_S1_C902_P99_N99_F#####.TIF`, **4256×4256 px**. Channel map from
  `RunSummary/Morphology_ChannelID_Dictionary.txt`: B=PanCK, G=G, Y=Membrane, R=CD45, U=DNA(DAPI).
- `CellStatsDir/FOV#####/` (182 dirs) — `CellBoundaries_F#####.csv`, `CellLabels_F#####.tif`,
  `CompartmentLabels_F#####.tif`, per-FOV cell stats.
- `CellStatsDir/CellComposite/`, `CellOverlay/`, `RnD/`,
  `Segmentation_2ad7413a-…_001/` (the segmentation set named in the metadata's
  `cellSegmentationSetId`).
- `AnalysisResults/6atoipi5h2/FOV#####/` — `…_complete_code_cell_target_call_coord.csv`, the
  per-FOV decoded target calls the flat files are derived from.
- `RunSummary/` — `…_ExptConfig.txt` (imaging params: `ImPixel_nm: 120.280945`, 8 cycles, 8 z-steps
  @ 0.8 µm, 4256² tiles), FOV locations, affine transform, spatial-BC metrics, QC/shading/distortion
  subdirs.
- `Logs/` — run and flow logs, scan params.

**Pixel size is 120.280945 nm** (`pixel_size_um = 0.120280945`); the sibling pipelines derive
`MM_PER_PX` from it. Global coordinates in the flat files are in this pixel unit.

## Running work on Sherlock

Per org policy, nothing heavy runs on the login node — no R, no Python, no image processing, no
whole-file `zcat` of the tx file. Use `sbatch` or `sh_dev`. This workspace already lives on
`$SCRATCH`, so job I/O belongs here (not `$HOME`); use `$L_SCRATCH` for high-IOPS temp and copy
results back before the job ends.

**R environment** — the group conda env, put on `PATH` rather than loaded via Lmod (there is no
`python` in it; R only):

```bash
export PATH="/oak/stanford/groups/bhowitt/conda_envs/workspace/bin:$PATH"
```

It carries Seurat v5, harmony, scPearsonPCA, gbspatial, clustree, dbscan, openxlsx, reticulate.
For anything else, `ml spider` first (`R/4.4.2` and `py-scanpy/1.10.2_py312` exist as modules).

**Partitions** — the lab owns `bhowitt` (2 nodes, 24 cores / 192 GB each, up to 7 d, use
`--qos=high_p`). It is small and often full; `owners` (up to 2 d) or `normal` are the fallbacks, and
`bigmem` for whole-object Seurat work. A full 126k-cell × 6.5k-feature Seurat object plus DecontX
needs well over the 8 GB/core default — request memory explicitly.

```bash
sbatch -p bhowitt --qos=high_p --time=12:00:00 --cpus-per-task=8 --mem=128GB job.sbatch
squeue --me
seff <jobid>
```

**`$SCRATCH` is purged after 90 days without a content write**, and this 39 GB export is not backed
up. Treat Box (`CosMX AIH/TA649/`, via `ml load system rclone`) or `$OAK` as the source of truth for
the raw data before relying on anything here long-term.
