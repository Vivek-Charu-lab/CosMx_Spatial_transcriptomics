## Spatial validation for the "stressed hepatocyte" call on clusters_round3 == G16:
## are G16 cells actually closer to immune infiltrate (T cells, plasma cells,
## macrophages, monocytes/neutrophils) than the healthy zonal hepatocyte
## clusters, consistent with AIH interface hepatitis rather than a random or
## purely technical subpopulation?
##
## Usage: Rscript src/g16_immune_proximity.R <seuratObj.rds> <celltype_table.tsv> <output_dir>

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(RANN)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
obj_path <- if (length(args) >= 1) args[1] else "data/after_despotx_seuratObj.rds"
celltype_path <- if (length(args) >= 2) args[2] else "celltype_table.tsv"
out_dir <- if (length(args) >= 3) args[3] else "results/stress_apoptosis"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hep_clusters <- c("G1", "G3", "G4", "G8", "G12", "G16")
immune_types <- c(
  "CD4_tem", "CD4_th17", "CD4_th2", "CD4_treg",
  "CD8_cytotoxic.1", "CD8_cytotoxic.2", "CD8_tem.1", "CD8_tem.2",
  "Macrophages_(Kupffer)", "Macrophages_(SAM)", "Monocytes.Neutrophils",
  "Plasma.cells.1_(Kappa)", "Plasma.cells.2_(Lambda)"
)

obj <- readRDS(obj_path)
sp <- Embeddings(obj, "spatial")
ct <- read.delim(celltype_path, stringsAsFactors = FALSE)

meta <- obj@meta.data %>%
  select(cell_id, clusters_round3, fov) %>%
  mutate(spatial_1 = sp[, 1], spatial_2 = sp[, 2]) %>%
  inner_join(ct, by = "cell_id")

immune_coords <- meta %>% filter(merged_celltypes %in% immune_types) %>%
  select(spatial_1, spatial_2) %>% as.matrix()
message("Immune cells for proximity query: ", nrow(immune_coords))

hep_meta <- meta %>% filter(clusters_round3 %in% hep_clusters)
hep_coords <- hep_meta %>% select(spatial_1, spatial_2) %>% as.matrix()

nn <- nn2(data = immune_coords, query = hep_coords, k = 1)
hep_meta$dist_to_nearest_immune <- as.vector(nn$nn.dists)

summary_tbl <- hep_meta %>%
  group_by(clusters_round3) %>%
  summarise(
    n_cells = n(),
    median_dist = median(dist_to_nearest_immune),
    mean_dist = mean(dist_to_nearest_immune),
    pct_within_50um = mean(dist_to_nearest_immune < 0.05) * 100,  # coords are in mm
    .groups = "drop"
  ) %>%
  arrange(median_dist)

write.table(summary_tbl, file.path(out_dir, "g16_immune_proximity_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

## Wilcoxon test: G16 vs. pooled healthy zonal hepatocytes (G1+G3+G4)
g16_d <- hep_meta$dist_to_nearest_immune[hep_meta$clusters_round3 == "G16"]
healthy_d <- hep_meta$dist_to_nearest_immune[hep_meta$clusters_round3 %in% c("G1", "G3", "G4")]
g8_d <- hep_meta$dist_to_nearest_immune[hep_meta$clusters_round3 == "G8"]
g12_d <- hep_meta$dist_to_nearest_immune[hep_meta$clusters_round3 == "G12"]

wtest_healthy <- wilcox.test(g16_d, healthy_d)
wtest_g8 <- wilcox.test(g16_d, g8_d)
wtest_g12 <- wilcox.test(g16_d, g12_d)

sink(file.path(out_dir, "g16_immune_proximity_stats.txt"))
cat("Distance to nearest immune cell (mm), by cluster:\n")
print(summary_tbl)
cat("\nWilcoxon rank-sum test, G16 vs healthy zonal hepatocytes (G1+G3+G4):\n")
print(wtest_healthy)
cat("\nWilcoxon rank-sum test, G16 vs G8 (Hep.IFNg):\n")
print(wtest_g8)
cat("\nWilcoxon rank-sum test, G16 vs G12 (Hep.Proliferating):\n")
print(wtest_g12)
sink()

hep_meta$clusters_round3 <- factor(hep_meta$clusters_round3, levels = summary_tbl$clusters_round3)
p <- ggplot(hep_meta, aes(x = clusters_round3, y = dist_to_nearest_immune, fill = clusters_round3)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  scale_y_continuous(limits = c(0, quantile(hep_meta$dist_to_nearest_immune, 0.98))) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Distance to nearest immune cell (T/plasma/macrophage/monocyte), by hepatocyte cluster",
       x = "clusters_round3", y = "distance (mm)")
ggsave(file.path(out_dir, "g16_immune_proximity.pdf"), p, width = 7, height = 5)

message("Done. Outputs in: ", normalizePath(out_dir))
print(summary_tbl)
