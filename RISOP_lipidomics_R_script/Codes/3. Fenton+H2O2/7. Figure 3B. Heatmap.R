rm(list = ls())
library(dplyr)
library(tidyverse)
library(stringr)
library(grid)
library(ComplexHeatmap)
library(circlize)
setwd(rstudioapi::getActiveProject())

# ============================================================
# 1. Load data # -----
# ============================================================

load(file = paste0("./Output/cmbd_lipids_Fenton.Rdata"))
load(file = paste0("./Output/cmbd_lipids_H2O2.Rdata"))

# ============================================================
# 2. Build OxPL sets and area tables # -----
# ============================================================

# 2.1 Filter OxPLs
cmbd_lipids_Fenton <- cmbd_lipids_Fenton %>%
  filter(grepl("<", Metabolite.curated, fixed = TRUE)) %>%
  select(Metabolite.curated)

cmbd_lipids_H2O2 <- cmbd_lipids_H2O2 %>%
  filter(grepl("<", Metabolite.curated, fixed = TRUE)) %>%
  select(Metabolite.curated)

# 2.2 Define shared and unique sets
fenton_set <- unique(cmbd_lipids_Fenton$Metabolite.curated)
h2o2_set   <- unique(cmbd_lipids_H2O2$Metabolite.curated)

Fenton_unique <- cmbd_lipids_Fenton %>% filter(!Metabolite.curated %in% h2o2_set)
H2O2_unique   <- cmbd_lipids_H2O2   %>% filter(!Metabolite.curated %in% fenton_set)
shared        <- cmbd_lipids_Fenton %>% filter(Metabolite.curated %in% h2o2_set)

# 2.3 Load peak area tables
Fenton_area <- read.csv('./Input/df_curated_Fenton.csv', stringsAsFactors = F, header = T) %>%
  filter(grepl("<", Metabolite.curated, fixed = TRUE)) %>%
  select(-Polarity, -Adduct.type, -Average.Rt.min., -Ontology)

H2O2_area <- read.csv('./Input/df_curated_H2O2.csv', stringsAsFactors = F, header = T) %>%
  filter(grepl("<", Metabolite.curated, fixed = TRUE)) %>%
  select(-Polarity, -Adduct.type, -Average.Rt.min., -Ontology)

# ============================================================
# 3. Normalize and annotate # -----
# ============================================================

# 3.1 Join area and normalize to sum
Fenton_unique <- Fenton_unique %>%
  left_join(Fenton_area, by = "Metabolite.curated") %>%
  mutate(across(-Metabolite.curated, ~ . / sum(., na.rm = TRUE))) %>%
  select(-matches("0h_"))

H2O2_unique <- H2O2_unique %>%
  left_join(H2O2_area, by = "Metabolite.curated") %>%
  mutate(across(-Metabolite.curated, ~ . / sum(., na.rm = TRUE))) %>%
  select(-matches("0h_"))

Fenton_shared <- shared %>%
  left_join(Fenton_area, by = "Metabolite.curated") %>%
  mutate(across(-Metabolite.curated, ~ . / sum(., na.rm = TRUE))) %>%
  select(-matches("0h_"))

H2O2_shared <- shared %>%
  left_join(H2O2_area, by = "Metabolite.curated") %>%
  mutate(across(-Metabolite.curated, ~ . / sum(., na.rm = TRUE))) %>%
  select(-matches("0h_"))

# 3.2 Count number of modifications per lipid
count_modifications <- function(mod) {
  if (is.na(mod)) return(NA_real_)
  parts <- str_trim(str_split(mod, ",")[[1]])
  total <- 0
  for (part in parts) {
    num   <- str_extract(part, "^\\d+")
    total <- total + ifelse(is.na(num), 1, as.numeric(num))
  }
  return(total)
}

# 3.3 Add modification annotation column
add_mod2 <- function(df) {
  df %>%
    mutate(Modification = str_extract(Metabolite.curated, "(?<=<)[^>]+(?=>)")) %>%
    mutate(
      mod_count     = sapply(Modification, count_modifications),
      Modification2 = case_when(
        is.na(Modification) ~ NA_character_,
        mod_count == 1      ~ paste0("<", Modification, ">"),
        mod_count > 1       ~ paste(mod_count, "modifications")
      )
    ) %>%
    select(-mod_count, -Modification) %>%
    mutate(Metabolite.curated = str_replace(Metabolite.curated,
                                            "^(\\w+)\\s+(.+?)\\s*(\\([a-z]\\))?$", "\\1(\\2)\\3"))
}

Fenton_unique <- add_mod2(Fenton_unique)
H2O2_unique   <- add_mod2(H2O2_unique)
Fenton_shared <- add_mod2(Fenton_shared)
H2O2_shared   <- add_mod2(H2O2_shared)

# ============================================================
# 4. Combine shared and unique sets # -----
# ============================================================

# 4.1 Merge shared Fenton and H2O2
Shared <- Fenton_shared %>%
  left_join(H2O2_shared, by = c("Metabolite.curated", "Modification2"))

shared_cols <- colnames(Shared)

# 4.2 Align unique sets to shared column structure
missing_in_Fenton <- setdiff(shared_cols, colnames(Fenton_unique))
Fenton_unique[missing_in_Fenton] <- NA
Fenton_unique <- Fenton_unique %>% select(all_of(shared_cols))

missing_in_H2O2 <- setdiff(shared_cols, colnames(H2O2_unique))
H2O2_unique[missing_in_H2O2] <- NA
H2O2_unique <- H2O2_unique %>% select(all_of(shared_cols))

# ============================================================
# 5. Average replicates and impute missing values # -----
# ============================================================

# 5.1 Average across time points per condition-replicate
conditions <- c("Fenton", "H2O2")
replicates <- c("1", "2", "3")

calc_mean <- function(df) {
  result <- df %>% select(Metabolite.curated, Modification2)
  for (cond in conditions) {
    for (rep in replicates) {
      cols <- grep(paste0("^", cond, "_.+_", rep, "$"), colnames(df), value = TRUE)
      if (length(cols) > 0) {
        result[[paste0(cond, "_", rep)]] <- rowMeans(df[, cols, drop = FALSE], na.rm = TRUE)
        result[[paste0(cond, "_", rep)]][is.nan(result[[paste0(cond, "_", rep)]])] <- NA
      } else {
        result[[paste0(cond, "_", rep)]] <- NA
      }
    }
  }
  return(result)
}

Shared_mean        <- calc_mean(Shared)
Fenton_unique_mean <- calc_mean(Fenton_unique)
H2O2_unique_mean   <- calc_mean(H2O2_unique)

# 5.2 Impute NAs as 5% of row minimum (unique sets only)
impute_min5pct <- function(df) {
  value_cols <- setdiff(colnames(df), c("Metabolite.curated", "Modification2"))
  df[value_cols] <- t(apply(df[value_cols], 1, function(x) {
    if (all(is.na(x))) return(x)
    row_min    <- min(x, na.rm = TRUE)
    x[is.na(x)] <- row_min * 0.05
    return(x)
  }))
  return(df)
}

Fenton_unique_mean <- impute_min5pct(Fenton_unique_mean)
H2O2_unique_mean   <- impute_min5pct(H2O2_unique_mean)

# ============================================================
# 6. Scale and prepare heatmap matrix # -----
# ============================================================

# 6.1 Bind all rows and z-score scale
value_cols <- c("Fenton_1", "Fenton_2", "Fenton_3", "H2O2_1", "H2O2_2", "H2O2_3")

df_full <- bind_rows(Shared_mean, Fenton_unique_mean, H2O2_unique_mean)

mat_full <- df_full %>%
  column_to_rownames("Metabolite.curated") %>%
  select(all_of(value_cols)) %>%
  as.matrix()

mat_scaled <- t(apply(mat_full, 1, function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x,   na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA, length(x)))
  (x - m) / s
}))
colnames(mat_scaled) <- value_cols

# 6.2 Define color scale
val_min <- min(mat_scaled, na.rm = TRUE)
val_max <- max(mat_scaled, na.rm = TRUE)
col_fun <- colorRamp2(c(val_min, 0, val_max), c("#2C48C7", "white", "#C91C4A"))

# ============================================================
# 7. Define heatmap aesthetics # -----
# ============================================================

# 7.1 Modification group levels, colors, and row title labels
mod_levels <- c("<oxo>", "<OH>", "<OOH>", "<COOH>", "2 modifications", "3 modifications")

group_colors <- c(
  "<oxo>"           = "#B4EEB4",
  "<OH>"            = "#FDB863",
  "<OOH>"           = "#7AC5CD",
  "<COOH>"          = "#EEB4B4",
  "2 modifications" = "#68838B",
  "3 modifications" = "#EED5B7"
)

row_title_labels <- c(
  "<oxo>"           = "Keto \u2212 <oxo>",
  "<OH>"            = "Hydroxy \u2212 <OH>",
  "<OOH>"           = "Hydroperoxy \u2212 <OOH>",
  "<COOH>"          = "Carboxy \u2212 <COOH>",
  "2 modifications" = "2 modifications",
  "3 modifications" = "3 modifications"
)

# 7.2 Column annotation (condition colors)
col_ha <- HeatmapAnnotation(
  Condition = anno_simple(
    c("Fenton", "Fenton", "Fenton", "H2O2", "H2O2", "H2O2"),
    col    = c("Fenton" = "#E05A51", "H2O2" = "#377EB8"),
    height = unit(3, "mm")
  ),
  show_annotation_name = FALSE,
  show_legend          = FALSE
)

# ============================================================
# 8. Build combined heatmap # -----
# ============================================================

# 8.1 Sort rows by modification and origin
pad <- 0.1

df_combined <- bind_rows(
  Shared_mean        %>% mutate(Origin = "Shared"),
  Fenton_unique_mean %>% mutate(Origin = "Fenton unique"),
  H2O2_unique_mean   %>% mutate(Origin = "H2O2 unique")
) %>%
  mutate(
    Modification2 = factor(Modification2, levels = mod_levels),
    Origin        = factor(Origin, levels = c("Shared", "Fenton unique", "H2O2 unique"))
  ) %>%
  arrange(Modification2, Origin)

mat_combined        <- mat_scaled[df_combined$Metabolite.curated, , drop = FALSE]
row_groups_combined <- df_combined$Modification2

# 8.2 Row annotation
row_ha_combined <- rowAnnotation(
  Modification = anno_simple(
    as.character(df_combined$Modification2),
    col   = group_colors,
    width = unit(3, "mm")
  ),
  show_annotation_name = FALSE
)

# 8.3 Cross-hatch cell function for not-detected entries
is_fenton_unique_row <- df_combined$Origin == "Fenton unique"
is_h2o2_unique_row   <- df_combined$Origin == "H2O2 unique"

cross_fun_combined <- function(j, i, x, y, width, height, fill) {
  need_cross <- (is_fenton_unique_row[i] && j %in% 4:6) ||
    (is_h2o2_unique_row[i]   && j %in% 1:3)
  if (need_cross) {
    grid.rect(x = x, y = y, width = width, height = height,
              gp = gpar(fill = "#D3D3D3", col = "white", lwd = 0.5))
    grid.lines(x = c(x - width*(0.5-pad), x + width*(0.5-pad)),
               y = c(y - height*(0.5-pad), y + height*(0.5-pad)),
               gp = gpar(col = "#5E5E5E", lwd = 0.5))
    grid.lines(x = c(x - width*(0.5-pad), x + width*(0.5-pad)),
               y = c(y + height*(0.5-pad), y - height*(0.5-pad)),
               gp = gpar(col = "#5E5E5E", lwd = 0.5))
  }
}

# 8.4 Build legends
nd_key_w <- unit(5, "mm")
nd_key_h <- unit(3, "mm")

zscore_legend_combined <- Legend(
  col_fun        = col_fun,
  title          = "Z score",
  title_gp       = gpar(fontsize = 8),
  labels_gp      = gpar(fontsize = 8),
  direction      = "horizontal",
  title_position = "lefttop",
  legend_width   = unit(24, "mm"),
  grid_height    = unit(3, "mm"),
  tick_length    = unit(0, "mm"),
  at             = c(val_min, 0, val_max),
  labels         = as.character(round(c(val_min, 0, val_max), 1))
)

not_detected_legend <- Legend(
  labels      = "Not detected",
  title       = "",
  title_gp    = gpar(fontsize = 0, col = "transparent"),
  labels_gp   = gpar(fontsize = 8),
  grid_width  = nd_key_w,
  grid_height = nd_key_h,
  graphics    = list(
    function(x, y, w, h) {
      grid.rect(x = x, y = y, width = nd_key_w, height = nd_key_h,
                gp = gpar(fill = "#D3D3D3", col = NA))
      half_w <- nd_key_w * (0.5 - pad)
      half_h <- nd_key_h * (0.5 - pad)
      grid.segments(x0 = x-half_w, y0 = y-half_h, x1 = x+half_w, y1 = y+half_h,
                    gp = gpar(col = "#5E5E5E", lwd = 0.8))
      grid.segments(x0 = x-half_w, y0 = y+half_h, x1 = x+half_w, y1 = y-half_h,
                    gp = gpar(col = "#5E5E5E", lwd = 0.8))
    }
  )
)

legend_packed <- packLegend(
  zscore_legend_combined,
  not_detected_legend,
  gap       = unit(6, "mm"),
  direction = "horizontal"
)

# 8.5 Assemble heatmap object
ht_combined <- Heatmap(
  mat_combined,
  name                = "Z score",
  col                 = col_fun,
  na_col              = "#D3D3D3",
  cell_fun            = cross_fun_combined,
  top_annotation      = col_ha,
  left_annotation     = row_ha_combined,
  cluster_columns     = FALSE,
  cluster_rows        = FALSE,
  row_split           = row_groups_combined,
  row_title           = row_title_labels[levels(row_groups_combined)],
  column_title_gp     = gpar(fontsize = 9),
  row_title_rot       = 0,
  row_title_gp        = gpar(fontsize = 9),
  show_row_names      = TRUE,
  show_column_names   = FALSE,
  row_names_gp        = gpar(fontsize = 5.5),
  rect_gp             = gpar(col = "white", lwd = 0.5),
  width               = unit(6, "cm"),
  show_heatmap_legend = FALSE
)

# ============================================================
# 9. Draw and export # -----
# ============================================================

draw_combined <- function() {
  draw(
    ht_combined,
    heatmap_legend_list    = list(legend_packed),
    heatmap_legend_side    = "bottom",
    show_annotation_legend = FALSE,
    merge_legend           = FALSE
  )
}

draw_combined()

pdf("./Figure/Heatmap_Combined_by_Modification.pdf", width = 830/96, height = 650/96)
draw_combined()
dev.off()