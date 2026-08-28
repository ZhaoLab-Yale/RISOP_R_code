rm(list = ls())
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(RColorBrewer)
setwd(rstudioapi::getActiveProject())

# ============================================================
# 1. Load data # -----
# ============================================================

load(file = "./Output/cmbd_lipids.Rdata")
load(file = "./Output/ML210_clusters.Rdata")
Sample_info <- read.csv('./Input/Sample_info.csv', stringsAsFactors = F, header = T)

# ============================================================
# 2. Pre-process lipid table # -----
# ============================================================

# 2.1 Filter to oxidized lipids and reformat names
cmbd_lipids <- cmbd_lipids %>%
  filter(grepl("<", Metabolite.curated, fixed = TRUE)) %>%
  mutate(Metabolite.curated = str_replace(Metabolite.curated,
                                          "^(\\w+)\\s+(.+?)\\s*(\\([a-z]\\))?$",
                                          "\\1(\\2)\\3"))

# 2.2 Build heatmap matrix
df_heatmap <- cmbd_lipids %>%
  select(-Ontology) %>%
  column_to_rownames(var = colnames(.)[1])

# ============================================================
# 3. Clean column names and define aesthetics # -----
# ============================================================

# 3.1 Reformat column names to "X min" style
colnames(df_heatmap) <- gsub("ML210_", "", colnames(df_heatmap))
colnames(df_heatmap) <- gsub("_",      " ", colnames(df_heatmap))
colnames(df_heatmap) <- gsub("(\\d+)(min)", "\\1 \\2", colnames(df_heatmap))

# 3.2 Timepoint and cluster color palettes
ann_colors <- list(
  Timepoint = c(
    "0 min"   = "#969696",
    "20 min"  = "#FFCFFD",
    "45 min"  = "#D296FF",
    "80 min"  = "#9591FF",
    "120 min" = "#3A5FCD",
    "270 min" = "#262175"
  )
)


# Cluster levels taken from the clustering result, no hard-coded k
cluster_levels <- sort(unique(as.character(ML210_clusters$Cluster)))
k              <- length(cluster_levels)

cluster_colors <- setNames(c("#A1863B", "#8B5742")[1:k], cluster_levels)

# 3.3 Map columns to timepoint factor
timepoints <- data.frame(
  Timepoint = factor(
    str_extract(colnames(df_heatmap), "\\d+ min"),
    levels = names(ann_colors$Timepoint)
  ),
  row.names = colnames(df_heatmap)
)

# ============================================================
# 4. Normalize cluster lipid names and build row order # -----
# ============================================================

# 4.1 Normalize lipid names to match heatmap row names
normalize_lipid_name <- function(name) {
  class_part  <- str_match(name, "^(\\w+) ")[, 2]
  middle_part <- str_match(name, "^\\w+ (.+?)( \\([a-z]\\))?$")[, 2]
  suffix_part <- str_match(name, " (\\([a-z]\\))$")[, 2]
  ifelse(!is.na(suffix_part),
         paste0(class_part, "(", middle_part, ")", suffix_part),
         paste0(class_part, "(", middle_part, ")"))
}

ML210_clusters_fixed <- ML210_clusters %>%
  mutate(Lipid   = normalize_lipid_name(as.character(Lipid)),
         Cluster = as.character(Cluster))

# 4.2 Join cluster assignments onto heatmap rows
row_cluster_df <- data.frame(Lipid = rownames(df_heatmap)) %>%
  left_join(ML210_clusters_fixed, by = "Lipid")

# 4.3 Sort rows by cluster then original order
row_cluster_df <- row_cluster_df %>%
  mutate(orig_order = row_number()) %>%
  arrange(factor(Cluster, levels = cluster_levels), orig_order)

df_heatmap_ordered <- df_heatmap[row_cluster_df$Lipid, ]
row_split          <- factor(row_cluster_df$Cluster, levels = cluster_levels)

# ============================================================
# 5. Order columns by timepoint and scale matrix # -----
# ============================================================

timepoints_ordered <- timepoints %>% arrange(Timepoint)
col_order          <- rownames(timepoints_ordered)
df_heatmap_final   <- df_heatmap_ordered[, col_order]
timepoints_final   <- timepoints[col_order, , drop = FALSE]

df_heatmap_final_scaled <- t(scale(t(df_heatmap_final)))
max_val <- max(abs(df_heatmap_final_scaled), na.rm = TRUE)

col_fun <- colorRamp2(c(-max_val, 0, max_val), c("blue", "white", "red"))

# ============================================================
# 6. Build heatmap annotations # -----
# ============================================================

# 6.1 Top annotation: timepoint label + color bar
tp_levels  <- levels(timepoints_final$Timepoint)
tp_col_vec <- as.character(timepoints_final$Timepoint)
col_split  <- factor(tp_col_vec, levels = tp_levels)

ha_final <- HeatmapAnnotation(
  Timepoint_label = anno_block(
    gp          = gpar(fill = NA, col = NA),
    labels      = tp_levels,
    labels_gp   = gpar(fontsize = 11, fontface = "bold", col = "black"),
    labels_just = "centre"
  ),
  Timepoint = timepoints_final$Timepoint,
  col = ann_colors,
  show_annotation_name = FALSE,
  show_legend          = FALSE,
  annotation_height    = unit(c(6, 5), "mm")
)

# 6.2 Left annotation: cluster color bar (width matches timepoint bar height)
ra_cluster <- rowAnnotation(
  Cluster = anno_block(
    gp        = gpar(fill = cluster_colors[levels(row_split)], col = NA),
    labels    = rep("", length(levels(row_split))),
    labels_gp = gpar(fontsize = 0)
  ),
  width = unit(5, "mm"),
  show_annotation_name = FALSE
)

# ============================================================
# 7. Build heatmap # -----
# ============================================================

row_titles <- paste("Cluster", levels(row_split))

ht_final <- Heatmap(df_heatmap_final_scaled,
                    name = "Z score",
                    col  = col_fun,
                    top_annotation    = ha_final,
                    left_annotation   = ra_cluster,
                    column_split      = col_split,
                    column_title      = NULL,
                    cluster_columns   = FALSE,
                    cluster_rows      = FALSE,
                    row_split         = row_split,
                    row_title         = row_titles,
                    row_title_gp      = gpar(fontsize = 11, fontface = "bold"),
                    row_gap           = unit(2, "mm"),
                    column_gap        = unit(0, "mm"),
                    border            = FALSE,
                    show_row_names    = TRUE,
                    show_column_names = FALSE,
                    row_names_gp      = gpar(fontsize = 9),
                    rect_gp           = gpar(col = "white", lwd = 0.5),
                    show_heatmap_legend = FALSE)

# ============================================================
# 8. Build Z-score legend # -----
# ============================================================

legend_zscore <- Legend(
  col_fun        = col_fun,
  title          = "Z score",
  direction      = "vertical",
  title_position = "topleft",
  title_gp       = gpar(fontsize = 9, fontface = "bold"),
  labels_gp      = gpar(fontsize = 9),
  grid_height    = unit(4,  "mm"),
  grid_width     = unit(4,  "mm"),
  legend_height  = unit(25, "mm")
)

# ============================================================
# 9. Draw and export # -----
# ============================================================

draw_heatmap <- function() {
  draw(ht_final,
       heatmap_legend_list    = list(legend_zscore),
       heatmap_legend_side    = "right",
       annotation_legend_side = "right",
       padding                = unit(c(2, 2, 2, 2), "mm"),
       align_heatmap_legend   = "heatmap_center"
  )
}

draw_heatmap()

pdf("./Figure/ML210_heatmap.pdf",
    width  = 700 / 96,
    height = 400 / 96)
draw_heatmap()
dev.off()