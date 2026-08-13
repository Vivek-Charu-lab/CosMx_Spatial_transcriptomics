## Stress / apoptosis module scoring on the post-DeSpotX TA649 Seurat object,
## evaluated per `clusters_round3`. Gene sets are literature-curated (MSigDB
## Hallmark Apoptosis/UPR/ROS, van den Brink/O'Flanagan/Denisenko core stress
## panel, AIH Fas/FasL apoptosis, hepatocyte identity) and pre-filtered to the
## TA649 6K Discovery panel.
##
## clusters_round3 is ~1 cluster per cell type (per celltype_table.tsv), not
## hepatocyte subclusters, so cross-cell-type score comparisons are not
## informative on their own: Hallmark Apoptosis includes ECM/stromal genes
## (BGN, DCN, LUM, TIMP1-3, MMP2) that read high in any fibroblast, and a
## generic "hepatocyte identity" score is trivially lowest in any non-hepatocyte
## cluster regardless of stress. The leaderboard below is therefore restricted
## to clusters whose dominant annotated cell type starts with "Hep".
##
## Usage: Rscript src/stress_apoptosis_scoring.R <path_to_after_despotx_seuratObj.rds> <path_to_celltype_table.tsv> <output_dir>

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
obj_path <- if (length(args) >= 1) args[1] else "data/after_despotx_seuratObj.rds"
celltype_path <- if (length(args) >= 2) args[2] else "celltype_table.tsv"
out_dir  <- if (length(args) >= 3) args[3] else "results/stress_apoptosis"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cluster_col <- "clusters_round3"

## ---- gene sets (literature) -------------------------------------------

gene_sets <- list(
  Stress_IEG_HSP = c(
    "FOS","FOSB","JUN","JUNB","JUND","ATF3","EGR1","NR4A1","IER3","BTG1",
    "BTG2","DUSP1","ZFP36","PPP1R15A","HSPA1A","HSPA1B","HSP90AA1",
    "HSP90AB1","HSPA8","HSPB1","DNAJB1"
  ),
  Apoptosis_Hallmark = c(
    "ADD1","AIFM3","ANKH","ANXA1","APP","ATF3","AVPR1A","BAX","BCAP31","BCL10",
    "BCL2L1","BCL2L10","BCL2L11","BCL2L2","BGN","BID","BIK","BIRC3","BMF","BMP2",
    "BNIP3L","BRCA1","BTG2","BTG3","CASP1","CASP2","CASP3","CASP4","CASP6","CASP7",
    "CASP8","CASP9","CAV1","CCNA1","CCND1","CCND2","CD14","CD2","CD38","CD44","CD69",
    "CDC25B","CDK2","CDKN1A","CDKN1B","CFLAR","CLU","CREBBP","CTH","CTNNB1","CYLD",
    "DAP","DAP3","DCN","DDIT3","DFFA","DIABLO","DNAJA1","DNAJC3","DNM1L","DPYD","EBP",
    "EGR3","EMP1","ENO2","ERBB2","ERBB3","EREG","ETF1","F2","F2R","FAS","FASLG","FDXR",
    "FEZ1","GADD45A","GADD45B","GCH1","GNA15","GPX1","GPX3","GPX4","GSN","GSR","GSTM1",
    "GUCY2D","H1F0","HGF","HMGB2","HMOX1","HSPB1","IER3","IFITM3","IFNB1","IFNGR1",
    "IGF2R","IGFBP6","IL18","IL1A","IL1B","IL6","IRF1","ISG20","JUN","KRT18","LEF1",
    "LGALS3","LMNA","LPPR4","LUM","MADD","MCL1","MGMT","MMP2","NEDD9","NEFH","PAK1",
    "PDCD4","PDGFRB","PEA15","PLAT","PLCB2","PMAIP1","PPP2R5B","PPP3R1","PPT1","PRF1",
    "PSEN1","PSEN2","PTK2","RARA","RELA","RETSAT","RHOB","RHOT2","RNASEL","ROCK1",
    "SAT1","SATB1","SC5DL","SLC20A1","SMAD7","SOD1","SOD2","SPTAN1","SQSTM1","TAP1",
    "TGFB2","TGFBR3","TIMP1","TIMP2","TIMP3","TNF","TNFRSF12A","TNFSF10","TOP2A","TSPO",
    "TXNIP","VDAC2","WEE1","XIAP"
  ),
  Apoptosis_Core_Mechanistic = c(
    "CASP3","CASP7","CASP8","CASP9","BAX","BAK1","BBC3","PMAIP1","BCL2L11","BID",
    "MCL1","BCL2L1","FAS","FASLG","TNFSF10","TNFRSF10B","XIAP","CFLAR","DFFA","DIABLO"
  ),
  Fas_FasL_AIH = c("FAS","FASLG","CFLAR","CASP8","CASP3"),
  UPR_Core = c(
    "ATF4","DDIT3","XBP1","ERN1","ATF6","EIF2AK3","DNAJC3","HERPUD1","HYOU1",
    "EDEM1","ERO1A"
  ),
  Oxidative_Stress_Core = c(
    "SOD1","SOD2","CAT","GPX1","GPX3","GPX4","GSR","GCLC","GCLM","TXNRD1",
    "PRDX1","NQO1","HMOX1","TXNIP"
  ),
  P53_Senescence_Core = c(
    "CDKN1A","GADD45A","GADD45B","MDM2","TP53","SESN1","TRIB3"
  ),
  Hepatocyte_Identity = c(
    "APOA1","APOA2","APOB","TTR","HNF4A","SERPINA1","CYP2E1"
  )
)

## ---- load & normalize ---------------------------------------------------

message("Loading object: ", obj_path)
obj <- readRDS(obj_path)
DefaultAssay(obj) <- "RNA"
stopifnot(cluster_col %in% colnames(obj@meta.data))

## ---- per-cluster cell-type identity (from independently annotated celltype_table) ----

celltype_df <- read.delim(celltype_path, stringsAsFactors = FALSE)
stopifnot(all(c("cell_id", "merged_celltypes") %in% colnames(celltype_df)))
cluster_composition <- obj@meta.data %>%
  select(cell_id, cluster = all_of(cluster_col)) %>%
  inner_join(celltype_df, by = "cell_id") %>%
  count(cluster, merged_celltypes) %>%
  group_by(cluster) %>%
  mutate(cluster_n = sum(n)) %>%
  ungroup()

cluster_celltype <- cluster_composition %>%
  group_by(cluster) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(cluster, dominant_celltype = merged_celltypes,
            dominant_pct = round(100 * n / cluster_n, 1))

write.table(cluster_composition %>% arrange(cluster, desc(n)),
            file.path(out_dir, "cluster_celltype_composition.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

if (!"data" %in% Layers(obj[["RNA"]]) ) {
  message("No normalized 'data' layer found — running NormalizeData (LogNormalize).")
  obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize",
                        scale.factor = 1e4, verbose = FALSE)
}

## ---- filter gene sets to genes actually on the panel/object -------------

present_genes <- rownames(obj)
gene_sets_filtered <- lapply(gene_sets, function(g) intersect(g, present_genes))
missing_report <- lapply(names(gene_sets), function(nm) {
  setdiff(gene_sets[[nm]], present_genes)
})
names(missing_report) <- names(gene_sets)

sink(file.path(out_dir, "gene_set_panel_coverage.txt"))
for (nm in names(gene_sets)) {
  cat(sprintf("%-28s present=%d/%d\n", nm, length(gene_sets_filtered[[nm]]),
              length(gene_sets[[nm]])))
  if (length(missing_report[[nm]]) > 0) {
    cat("  missing: ", paste(missing_report[[nm]], collapse = ", "), "\n")
  }
}
sink()

## ---- AddModuleScore -------------------------------------------------------

score_names <- names(gene_sets_filtered)
obj <- AddModuleScore(
  obj,
  features = gene_sets_filtered,
  name = paste0(score_names, "_"),
  seed = 1,
  nbin = 24
)

## AddModuleScore appends the running list-index (not always "1") to each
## name[i] prefix; recover mapping and rename to clean score names.
score_cols <- paste0(score_names, "_", seq_along(score_names))
stopifnot(all(score_cols %in% colnames(obj@meta.data)))
for (i in seq_along(score_names)) {
  obj@meta.data[[score_names[i]]] <- obj@meta.data[[score_cols[i]]]
}
obj@meta.data[score_cols] <- NULL

## ---- per-cluster summary ---------------------------------------------------

meta <- obj@meta.data
meta$cluster <- meta[[cluster_col]]

summary_tbl <- meta %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    median_nCount = median(nCount_RNA),
    median_nFeature = median(nFeature_RNA),
    across(all_of(score_names), list(mean = ~mean(.x), median = ~median(.x)),
           .names = "{.col}__{.fn}"),
    .groups = "drop"
  ) %>%
  arrange(desc(Stress_IEG_HSP__mean))

summary_tbl <- summary_tbl %>%
  left_join(cluster_celltype, by = "cluster") %>%
  relocate(dominant_celltype, dominant_pct, .after = cluster)

write.table(summary_tbl, file.path(out_dir, "module_scores_by_cluster.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

## Rank clusters jointly by stress + apoptosis + identity-loss. This is only
## meaningful WITHIN hepatocyte clusters (dominant_celltype starts with "Hep"
## at >=90% purity) -- across cell types it just re-derives "is this a
## hepatocyte", not "is this hepatocyte stressed".
joint_rank_all <- meta %>%
  group_by(cluster) %>%
  summarise(
    mean_stress = mean(Stress_IEG_HSP),
    mean_apoptosis_hallmark = mean(Apoptosis_Hallmark),
    mean_apoptosis_core = mean(Apoptosis_Core_Mechanistic),
    mean_fas_fasl = mean(Fas_FasL_AIH),
    mean_upr = mean(UPR_Core),
    mean_oxidative = mean(Oxidative_Stress_Core),
    mean_p53 = mean(P53_Senescence_Core),
    mean_identity = mean(Hepatocyte_Identity),
    median_nCount = median(nCount_RNA),
    median_nFeature = median(nFeature_RNA),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  left_join(cluster_celltype, by = "cluster") %>%
  relocate(dominant_celltype, dominant_pct, .after = cluster)

write.table(joint_rank_all %>% arrange(desc(mean_stress)),
            file.path(out_dir, "cluster_scores_all_celltypes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

hepatocyte_clusters <- joint_rank_all %>%
  filter(startsWith(dominant_celltype, "Hep"), dominant_pct >= 90) %>%
  pull(cluster)
message("Hepatocyte clusters used for leaderboard: ", paste(hepatocyte_clusters, collapse = ", "))

joint_rank <- joint_rank_all %>%
  filter(cluster %in% hepatocyte_clusters) %>%
  mutate(
    stress_rank = rank(-mean_stress),
    apoptosis_rank = rank(-mean_apoptosis_hallmark),
    apoptosis_core_rank = rank(-mean_apoptosis_core),
    identity_loss_rank = rank(mean_identity),
    joint_rank_score = stress_rank + apoptosis_rank + apoptosis_core_rank
  ) %>%
  arrange(joint_rank_score)

write.table(joint_rank, file.path(out_dir, "hepatocyte_cluster_ranking.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

## ---- plots ---------------------------------------------------------------

## (1) all clusters, labeled by dominant cell type, for transparency about
## what clusters_round3 actually contains.
celltype_label <- setNames(
  paste0(cluster_celltype$cluster, ": ", cluster_celltype$dominant_celltype,
         " (", cluster_celltype$dominant_pct, "%)"),
  cluster_celltype$cluster
)
meta_all <- meta
meta_all$cluster_label <- factor(celltype_label[as.character(meta_all$cluster)],
                                  levels = celltype_label[levels(factor(meta_all$cluster))])

pdf(file.path(out_dir, "module_score_violins_all_celltypes.pdf"), width = 13, height = 7)
for (sn in score_names) {
  p <- ggplot(meta_all, aes(x = cluster_label, y = .data[[sn]], fill = cluster_label)) +
    geom_violin(scale = "width", trim = TRUE) +
    stat_summary(fun = median, geom = "point", size = 1, color = "black") +
    theme_minimal() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0(sn, " (all clusters_round3, by dominant cell type)"),
         x = NULL, y = "AddModuleScore")
  print(p)
}
dev.off()

## (2) hepatocyte-only leaderboard, ranked by joint_rank_score.
meta_hep <- meta %>% filter(cluster %in% hepatocyte_clusters)
meta_hep$cluster <- factor(meta_hep$cluster, levels = joint_rank$cluster)

pdf(file.path(out_dir, "module_score_violins_hepatocytes_only.pdf"), width = 9, height = 6)
for (sn in score_names) {
  p <- ggplot(meta_hep, aes(x = cluster, y = .data[[sn]], fill = cluster)) +
    geom_violin(scale = "width", trim = TRUE) +
    stat_summary(fun = median, geom = "point", size = 1, color = "black") +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = paste0(sn, " (hepatocyte clusters only, ranked by joint stress+apoptosis score)"),
         x = "clusters_round3", y = "AddModuleScore")
  print(p)
}
dev.off()

## Heatmap of mean module score (z-scored across the hepatocyte clusters only).
score_mat <- as.matrix(joint_rank[, c("mean_stress","mean_apoptosis_hallmark",
                                       "mean_apoptosis_core","mean_fas_fasl",
                                       "mean_upr","mean_oxidative","mean_p53",
                                       "mean_identity")])
rownames(score_mat) <- paste0(joint_rank$cluster, " (", joint_rank$dominant_celltype, ")")
z_mat <- scale(score_mat)

pdf(file.path(out_dir, "module_score_heatmap_hepatocytes_only.pdf"), width = 9, height = 6)
heatmap(z_mat, Colv = NA, scale = "none", margins = c(10, 14),
        main = "Mean module score (z-scored across hepatocyte clusters)")
dev.off()

message("Done. Outputs written to: ", normalizePath(out_dir))
cat("\n--- Hepatocyte-cluster leaderboard (stress + apoptosis) ---\n")
print(joint_rank %>% select(cluster, dominant_celltype, dominant_pct, mean_stress,
                             mean_apoptosis_hallmark, mean_apoptosis_core, mean_identity,
                             joint_rank_score))
cat("\n--- All clusters, for context (cell type composition) ---\n")
print(joint_rank_all %>% select(cluster, dominant_celltype, dominant_pct, n_cells) %>%
        arrange(cluster))
