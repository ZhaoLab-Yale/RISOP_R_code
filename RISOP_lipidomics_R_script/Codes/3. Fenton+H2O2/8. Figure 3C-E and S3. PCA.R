rm(list = ls())
library(tidyverse)
library(openxlsx)
library(ggrepel)
library(svglite)
library(cowplot)
setwd(rstudioapi::getActiveProject())

# ============================================================
# 1. Load data # -----
# ============================================================

load(file = paste0("./Output/cmbd_lipids_Fenton.Rdata"))
load(file = paste0("./Output/cmbd_lipids_H2O2.Rdata"))
Sample_info <- read.csv('./Input/Sample_info.csv', stringsAsFactors = F, skip = 0, header = T)

# ============================================================
# 2. Shared theme, colors, shapes, and labels # -----
# ============================================================

# 2.1 ggplot theme
JCtheme_PCA <- function() {
  theme_bw() +
    theme(
      panel.background      = element_blank(),
      panel.grid.major      = element_blank(),
      panel.grid.minor      = element_blank(),
      panel.border          = element_rect(color = "black", fill = NA, size = 1),
      axis.text.x           = element_text(angle = 0, hjust = 0.5, vjust = 1,
                                           size = 16, color = "black",
                                           face = "plain", family = "Arial"),
      axis.text.y           = element_text(hjust = 1, vjust = 0.5,
                                           size = 16, color = "black",
                                           face = "plain", family = "Arial"),
      axis.title.x          = element_text(angle = 0, hjust = 0.5, vjust = -1.5,
                                           face = "bold", family = "Arial",
                                           colour = "black", size = 18),
      axis.title.y          = element_text(angle = 90, hjust = 0.5, vjust = 3,
                                           face = "bold", family = "Arial",
                                           colour = "black", size = 18),
      legend.text           = element_text(face = "plain", family = "Arial",
                                           colour = "black", size = 14.5),
      legend.title          = element_text(face = "bold", family = "Arial",
                                           colour = "black", size = 15),
      legend.title.align    = 0.27,
      legend.position       = "right",
      legend.justification  = c(0, 0),
      legend.box.just       = "bottom",
      legend.key.width      = unit(0.6, "cm"),
      legend.key.height     = unit(0.6, "cm"),
      legend.box.spacing    = unit(0.35, "cm"),
      plot.margin           = margin(t = 5, b = 5, l = 5, r = 5),
      aspect.ratio          = 1
    )
}

# 2.2 Combined (Fenton + H2O2) colors, shapes, and labels
color_values_combined <- c(
  "Fenton_0 h"    = "#969696",
  "Fenton_1 min"  = "#FCBBA1",
  "Fenton_20 min" = "#FC9272",
  "Fenton_4 h"    = "#FB6A4A",
  "Fenton_24 h"   = "#DE2D26",
  "Fenton_72 h"   = "#A50F15",
  "H2O2_0 h"      = "#969696",
  "H2O2_1 min"    = "#C6DBEF",
  "H2O2_20 min"   = "#9ECAE1",
  "H2O2_4 h"      = "#6BAED6",
  "H2O2_24 h"     = "#3182BD",
  "H2O2_72 h"     = "#08519C"
)

shape_values_combined <- c(
  "Fenton_0 h"    = 16, "Fenton_1 min"  = 16, "Fenton_20 min" = 16,
  "Fenton_4 h"    = 16, "Fenton_24 h"   = 16, "Fenton_72 h"   = 16,
  "H2O2_0 h"      = 17, "H2O2_1 min"    = 17, "H2O2_20 min"   = 17,
  "H2O2_4 h"      = 17, "H2O2_24 h"     = 17, "H2O2_72 h"     = 17
)

custom_labels_combined <- c(
  "Fenton_0 h"    = "Fenton 0 h",
  "Fenton_1 min"  = "Fenton 1 min",
  "Fenton_20 min" = "Fenton 20 min",
  "Fenton_4 h"    = "Fenton 4 h",
  "Fenton_24 h"   = "Fenton 24 h",
  "Fenton_72 h"   = "Fenton 72 h",
  "H2O2_0 h"      = expression(plain(H[2]*O[2]*"-only 0 h")),
  "H2O2_1 min"    = expression(plain(H[2]*O[2]*"-only 1 min")),
  "H2O2_20 min"   = expression(plain(H[2]*O[2]*"-only 20 min")),
  "H2O2_4 h"      = expression(plain(H[2]*O[2]*"-only 4 h")),
  "H2O2_24 h"     = expression(plain(H[2]*O[2]*"-only 24 h")),
  "H2O2_72 h"     = expression(plain(H[2]*O[2]*"-only 72 h"))
)

# 2.3 H2O2-only colors and labels
color_values_H2O2 <- c(
  "0 h"    = "#969696",
  "1 min"  = "#C6DBEF",
  "20 min" = "#9ECAE1",
  "4 h"    = "#6BAED6",
  "24 h"   = "#3182BD",
  "72 h"   = "#08519C"
)

custom_labels_H2O2 <- c(
  "0 h"    = expression(plain(H[2]*O[2]*"-only 0 h")),
  "1 min"  = expression(plain(H[2]*O[2]*"-only 1 min")),
  "20 min" = expression(plain(H[2]*O[2]*"-only 20 min")),
  "4 h"    = expression(plain(H[2]*O[2]*"-only 4 h")),
  "24 h"   = expression(plain(H[2]*O[2]*"-only 24 h")),
  "72 h"   = expression(plain(H[2]*O[2]*"-only 72 h"))
)

# 2.4 Fenton-only colors and labels
color_values_Fenton <- c(
  "0 h"    = "#969696",
  "1 min"  = "#FCBBA1",
  "20 min" = "#FC9272",
  "4 h"    = "#FB6A4A",
  "24 h"   = "#DE2D26",
  "72 h"   = "#A50F15"
)

custom_labels_Fenton <- c(
  "0 h"    = "Fenton 0 h",
  "1 min"  = "Fenton 1 min",
  "20 min" = "Fenton 20 min",
  "4 h"    = "Fenton 4 h",
  "24 h"   = "Fenton 24 h",
  "72 h"   = "Fenton 72 h"
)

# ============================================================
# 3. PCA functions # -----
# ============================================================

# 3.1 Run PCA for a single dataset (colored by Time)
run_pca_single <- function(cmbd_lipids, Sample_info) {
  df <- cmbd_lipids %>%
    ungroup() %>%
    select(-Ontology) %>%
    column_to_rownames(var = colnames(.)[1]) %>%
    t() %>%
    as.data.frame() %>%
    log2()
  
  df_Group <- df %>%
    rownames_to_column(var = "Sample_name") %>%
    left_join(Sample_info %>% select(Sample_name, Time), by = "Sample_name") %>%
    column_to_rownames(var = "Sample_name")
  
  df_pca   <- prcomp(df, center = TRUE, scale. = TRUE)
  eig      <- (df_pca$sdev)^2
  variance <- eig * 100 / sum(eig)
  
  df_pcs <- data.frame(df_pca$x,
                       Group       = df_Group$Time,
                       Sample_name = rownames(df_Group))
  
  df_pcs$Group <- factor(df_pcs$Group,
                         levels = c("0 h", "1 min", "20 min", "4 h", "24 h", "72 h"))
  
  list(df_pcs   = df_pcs,
       variance = variance)
}

# 3.2 Run PCA for combined dataset (colored by Group_Time)
run_pca_combined <- function(df_features, Sample_info) {
  df <- df_features %>%
    column_to_rownames(var = colnames(.)[1]) %>%
    t() %>%
    as.data.frame() %>%
    log2()
  
  df_Group <- df %>%
    rownames_to_column(var = "Sample_name") %>%
    left_join(Sample_info %>% select(Sample_name, Group, Time),
              by = "Sample_name") %>%
    column_to_rownames(var = "Sample_name")
  
  df_pca   <- prcomp(df, center = TRUE, scale. = TRUE)
  eig      <- (df_pca$sdev)^2
  variance <- eig * 100 / sum(eig)
  
  df_pcs <- data.frame(df_pca$x,
                       Group       = df_Group$Group,
                       Time        = df_Group$Time,
                       Sample_name = rownames(df_Group))
  
  df_pcs$Time <- factor(df_pcs$Time,
                        levels = c("0 h", "1 min", "20 min", "4 h", "24 h", "72 h"))
  
  df_pcs <- df_pcs %>%
    mutate(Group_Time = factor(paste(Group, Time, sep = "_"),
                               levels = names(color_values_combined)))
  
  list(df_pcs   = df_pcs,
       variance = variance)
}

# 3.3 Build PCA plot for a single dataset
make_pca_plot_single <- function(pca_result, color_values, color_labels, shape_val) {
  df_pcs     <- pca_result$df_pcs
  variance   <- pca_result$variance
  xaxistitle <- paste0("PC1: ", round(variance[1], 2), "%")
  yaxistitle <- paste0("PC2: ", round(variance[2], 2), "%")
  
  ggplot(df_pcs, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
    geom_hline(aes(yintercept = 0), color = "black", lty = "dashed") +
    geom_vline(aes(xintercept = 0), color = "black", lty = "dashed") +
    geom_point(size = 4) +
    scale_color_manual(name   = "Conditions",
                       values = color_values,
                       labels = color_labels) +
    scale_shape_manual(name   = "Conditions",
                       values = rep(shape_val, 6),
                       labels = color_labels) +
    labs(x = xaxistitle, y = yaxistitle) +
    JCtheme_PCA() +
    guides(
      color = guide_legend(override.aes = list(size = 3.5)),
      shape = guide_legend(override.aes = list(size = 3.5))
    )
}

# 3.4 Build PCA plot for combined dataset
make_pca_plot_combined <- function(pca_result) {
  df_pcs     <- pca_result$df_pcs
  variance   <- pca_result$variance
  xaxistitle <- paste0("PC1: ", round(variance[1], 2), "%")
  yaxistitle <- paste0("PC2: ", round(variance[2], 2), "%")
  
  ggplot(df_pcs, aes(x = PC1, y = PC2,
                     color = Group_Time, shape = Group_Time)) +
    geom_hline(aes(yintercept = 0), color = "black", lty = "dashed") +
    geom_vline(aes(xintercept = 0), color = "black", lty = "dashed") +
    geom_point(size = 4) +
    labs(x = xaxistitle, y = yaxistitle) +
    scale_color_manual(name   = "Conditions",
                       values = color_values_combined,
                       labels = custom_labels_combined) +
    scale_shape_manual(name   = "Conditions",
                       values = shape_values_combined,
                       labels = custom_labels_combined) +
    JCtheme_PCA() +
    guides(
      color = guide_legend(override.aes = list(size = 3.5)),
      shape = guide_legend(override.aes = list(size = 3.5))
    )
}

# ============================================================
# 4. Prepare combined (Fenton + H2O2) dataset # -----
# ============================================================

shared_metabolites <- intersect(cmbd_lipids_Fenton$Metabolite.curated,
                                cmbd_lipids_H2O2$Metabolite.curated)

cmbd_lipids_Fenton_shared <- cmbd_lipids_Fenton %>%
  filter(Metabolite.curated %in% shared_metabolites) %>%
  select(-Ontology)

cmbd_lipids_H2O2_shared <- cmbd_lipids_H2O2 %>%
  filter(Metabolite.curated %in% shared_metabolites) %>%
  select(-Ontology)

df_curated <- cmbd_lipids_Fenton_shared %>%
  left_join(cmbd_lipids_H2O2_shared, by = "Metabolite.curated")

# ============================================================
# 5. PCA input feature counts # -----
# ============================================================

count_features <- function(df, label) {
  n_lipids  <- nrow(df)
  meta_cols <- c("Metabolite.curated", "Ontology")
  n_samples <- ncol(df) - sum(colnames(df) %in% meta_cols)
  message(sprintf("[%-35s]  lipids = %d,  samples = %d", label, n_lipids, n_samples))
}

message("\n========== PCA Input Feature Counts ==========")

count_features(
  df_curated %>% filter( grepl("<", Metabolite.curated, fixed = TRUE)),
  "Combined Oxidized"
)
count_features(
  df_curated %>% filter(!grepl("<", Metabolite.curated, fixed = TRUE)),
  "Combined Non-oxidized"
)
count_features(
  cmbd_lipids_H2O2 %>% filter( grepl("<", Metabolite.curated, fixed = TRUE)),
  "H2O2 Oxidized"
)
count_features(
  cmbd_lipids_H2O2 %>% filter(!grepl("<", Metabolite.curated, fixed = TRUE)),
  "H2O2 Non-oxidized"
)
count_features(
  cmbd_lipids_Fenton %>% filter( grepl("<", Metabolite.curated, fixed = TRUE)),
  "Fenton Oxidized"
)
count_features(
  cmbd_lipids_Fenton %>% filter(!grepl("<", Metabolite.curated, fixed = TRUE)),
  "Fenton Non-oxidized"
)

message("===============================================\n")

# ============================================================
# 6. Run all 6 PCAs # -----
# ============================================================

pca_combined_Ox    <- run_pca_combined(df_curated %>% filter( grepl("<", Metabolite.curated, fixed = TRUE)), Sample_info)
pca_combined_nonOx <- run_pca_combined(df_curated %>% filter(!grepl("<", Metabolite.curated, fixed = TRUE)), Sample_info)

pca_H2O2_Ox    <- run_pca_single(cmbd_lipids_H2O2 %>% filter( grepl("<", Metabolite.curated, fixed = TRUE)), Sample_info)
pca_H2O2_nonOx <- run_pca_single(cmbd_lipids_H2O2 %>% filter(!grepl("<", Metabolite.curated, fixed = TRUE)), Sample_info)

pca_Fenton_Ox    <- run_pca_single(cmbd_lipids_Fenton %>% filter( grepl("<", Metabolite.curated, fixed = TRUE)), Sample_info)
pca_Fenton_nonOx <- run_pca_single(cmbd_lipids_Fenton %>% filter(!grepl("<", Metabolite.curated, fixed = TRUE)), Sample_info)

# ============================================================
# 7. Build all 6 PCA plots # -----
# ============================================================

plot_combined_Ox    <- make_pca_plot_combined(pca_combined_Ox)
plot_combined_nonOx <- make_pca_plot_combined(pca_combined_nonOx)

plot_H2O2_Ox    <- make_pca_plot_single(pca_H2O2_Ox,    color_values_H2O2,   custom_labels_H2O2,   shape_val = 17)
plot_H2O2_nonOx <- make_pca_plot_single(pca_H2O2_nonOx, color_values_H2O2,   custom_labels_H2O2,   shape_val = 17)

plot_Fenton_Ox    <- make_pca_plot_single(pca_Fenton_Ox,    color_values_Fenton, custom_labels_Fenton, shape_val = 16)
plot_Fenton_nonOx <- make_pca_plot_single(pca_Fenton_nonOx, color_values_Fenton, custom_labels_Fenton, shape_val = 16)

# ============================================================
# 8. Align and export all 6 plots # -----
# ============================================================

# 8.1 Align plots
aligned <- align_plots(
  plot_combined_Ox,
  plot_combined_nonOx,
  plot_H2O2_Ox,
  plot_H2O2_nonOx,
  plot_Fenton_Ox,
  plot_Fenton_nonOx,
  align = "hv",
  axis  = "tblr"
)

# 8.2 Export as SVG
ggsave("./Figure/Ox_PCA_Fenton_H2O2.svg",
       plot   = aligned[[1]],
       width  = 610 / 96, height = 440 / 96, units = "in",
       device = svglite::svglite, bg = "transparent")

ggsave("./Figure/nonOx_PCA_Fenton_H2O2.svg",
       plot   = aligned[[2]],
       width  = 610 / 96, height = 440 / 96, units = "in",
       device = svglite::svglite, bg = "transparent")

ggsave("./Figure/Ox_H2O2_PCA.svg",
       plot   = aligned[[3]],
       width  = 610 / 96, height = 440 / 96, units = "in",
       device = svglite::svglite, bg = "transparent")

ggsave("./Figure/nonOx_H2O2_PCA.svg",
       plot   = aligned[[4]],
       width  = 610 / 96, height = 440 / 96, units = "in",
       device = svglite::svglite, bg = "transparent")

ggsave("./Figure/Ox_Fenton_PCA.svg",
       plot   = aligned[[5]],
       width  = 610 / 96, height = 440 / 96, units = "in",
       device = svglite::svglite, bg = "transparent")

ggsave("./Figure/nonOx_Fenton_PCA.svg",
       plot   = aligned[[6]],
       width  = 610 / 96, height = 440 / 96, units = "in",
       device = svglite::svglite, bg = "transparent")