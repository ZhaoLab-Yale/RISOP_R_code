rm(list = ls())
library(tidyverse)
library(ggrepel)
setwd(rstudioapi::getActiveProject())

# ============================================================
# 1. Load data # -----
# ============================================================

load(file = "./Output/cmbd_lipids.Rdata")
Sample_info <- read.csv('./Input/Sample_info.csv', stringsAsFactors = F, skip = 0, header = T)

# ============================================================
# 2. Shared theme, colors, and time-point levels # -----
# ============================================================

# 2.1 Time-point levels and color gradient
time_levels   <- c("0 min", "20 min", "45 min", "80 min", "120 min", "270 min")
teal_gradient <- c("#969696", "#FFCFFD", "#D296FF", "#9591FF", "#3A5FCD", "#262175")

# 2.2 ggplot theme
JCtheme_PCA <- function() {
  theme_bw() +
    theme(
      panel.background     = element_blank(),
      panel.grid.major     = element_blank(),
      panel.grid.minor     = element_blank(),
      panel.border         = element_rect(color = "black", fill = NA, size = 1),
      axis.text.x          = element_text(angle = 0, hjust = 0.5, vjust = 1,
                                          size = 14, color = "black",
                                          face = "plain", family = "Arial"),
      axis.text.y          = element_text(hjust = 1, vjust = 0.5,
                                          size = 14, color = "black",
                                          face = "plain", family = "Arial"),
      axis.title.x         = element_text(angle = 0, hjust = 0.5, vjust = -1.5,
                                          face = "bold", family = "Arial",
                                          colour = "black", size = 16),
      axis.title.y         = element_text(angle = 90, hjust = 0.5, vjust = 3,
                                          face = "bold", family = "Arial",
                                          colour = "black", size = 16),
      legend.text          = element_text(face = "plain", family = "Arial",
                                          colour = "black", size = 14),
      legend.title         = element_text(face = "bold", family = "Arial",
                                          colour = "black", size = 15),
      legend.title.align   = 0.27,
      legend.justification = c(0, -0.12),
      legend.box.spacing   = unit(0.35, "cm"),
      plot.title           = element_text(size = 18, face = "bold",
                                          family = "Arial", hjust = 0.5),
      plot.margin          = margin(t = 20, b = 15, l = 20, r = 25),
      aspect.ratio         = 1
    )
}

# ============================================================
# 3. PCA functions # -----
# ============================================================

# 3.1 Run PCA
run_pca_ml210 <- function(cmbd_lipids, Sample_info) {
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
  
  df_pcs$Group <- factor(df_pcs$Group, levels = time_levels)
  
  xaxistitle <- paste0("PC1: ", round(variance[1], 2), "%")
  yaxistitle <- paste0("PC2: ", round(variance[2], 2), "%")
  
  list(df_pcs     = df_pcs,
       xaxistitle = xaxistitle,
       yaxistitle = yaxistitle)
}

# 3.2 Build PCA plot
make_pca_plot_ml210 <- function(pca_result) {
  df_pcs <- pca_result$df_pcs
  
  ggplot(df_pcs, aes(x = PC1, y = PC2, color = Group, shape = Group, fill = Group)) +
    geom_hline(aes(yintercept = 0), color = "black", lty = "dashed") +
    geom_vline(aes(xintercept = 0), color = "black", lty = "dashed") +
    geom_point(size = 4) +
    scale_color_manual(name   = "Time points",
                       values = teal_gradient,
                       labels = time_levels,
                       drop   = FALSE) +
    scale_shape_manual(name   = "Time points",
                       values = rep(23, 6),
                       labels = time_levels,
                       drop   = FALSE) +
    scale_fill_manual(name   = "Time points",
                      values = teal_gradient,
                      labels = time_levels,
                      drop   = FALSE) +
    labs(x = pca_result$xaxistitle, y = pca_result$yaxistitle) +
    JCtheme_PCA() +
    guides(
      color = guide_legend(override.aes = list(size = 3.5)),
      shape = guide_legend(override.aes = list(size = 3.5)),
      fill  = guide_legend(override.aes = list(size = 3.5))
    )
}

# ============================================================
# 4. Non-oxidized lipids PCA # -----
# ============================================================

cmbd_lipids_nonOx <- cmbd_lipids %>%
  filter(!grepl("<", Metabolite.curated, fixed = TRUE))

message("Non-oxidized lipids: ", nrow(cmbd_lipids_nonOx))

pca_nonOx  <- run_pca_ml210(cmbd_lipids_nonOx, Sample_info)
plot_nonOx <- make_pca_plot_ml210(pca_nonOx)

print(plot_nonOx)

ggsave("./Figure/NonOx_ML210_PCA.svg",
       plot   = plot_nonOx,
       width  = 508 / 96, height = 410 / 96, units = "in",
       device = svglite::svglite,
       bg     = "transparent")

# ============================================================
# 5. Oxidized lipids PCA # -----
# ============================================================

cmbd_lipids_Ox <- cmbd_lipids %>%
  filter(grepl("<", Metabolite.curated, fixed = TRUE))

message("Oxidized lipids: ", nrow(cmbd_lipids_Ox))

pca_Ox  <- run_pca_ml210(cmbd_lipids_Ox, Sample_info)
plot_Ox <- make_pca_plot_ml210(pca_Ox)

print(plot_Ox)

ggsave("./Figure/Ox_ML210_PCA.svg",
       plot   = plot_Ox,
       width  = 508 / 96, height = 410 / 96, units = "in",
       device = svglite::svglite,
       bg     = "transparent")