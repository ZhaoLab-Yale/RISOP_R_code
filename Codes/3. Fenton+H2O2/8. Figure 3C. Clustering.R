rm(list = ls())
library(tidyverse)
library(gridExtra)
library(grid)
library(ggpubr)
library(svglite)
library(fpc)
setwd(rstudioapi::getActiveProject())

# ============================================================
# 0. Global settings # -----
# ============================================================

K_SCAN     <- 2:9    # candidate k values to assess for stability
BOOT_B     <- 1000   # bootstrap iterations
KM_NSTART  <- 25     # random starts for K-means; also passed to clusterboot
SEED       <- 12345    # single seed used for both clusterboot and K-means
K_CLUSTERS <- 4      # k used for the final figure and saved clusters

# ============================================================
# 1. Load data # -----
# ============================================================

load(file = paste0("./Output/cmbd_lipids_Fenton.Rdata"))
load(file = paste0("./Output/cmbd_lipids_H2O2.Rdata"))
Sample_info <- read.csv('./Input/Sample_info.csv', stringsAsFactors = F, skip = 0, header = T)

# ============================================================
# 2. Filter to shared lipids and stack datasets # -----
# ============================================================

# 2.1 Identify shared metabolites
shared_metabolites <- intersect(cmbd_lipids_Fenton$Metabolite.curated,
                                cmbd_lipids_H2O2$Metabolite.curated)

cmbd_lipids_Fenton <- cmbd_lipids_Fenton %>%
  filter(Metabolite.curated %in% shared_metabolites)

cmbd_lipids_H2O2 <- cmbd_lipids_H2O2 %>%
  filter(Metabolite.curated %in% shared_metabolites)

# 2.2 Append treatment label and strip column prefixes before stacking
df_Fenton <- cmbd_lipids_Fenton %>%
  mutate(Metabolite.curated = paste0(Metabolite.curated, " (Fenton)")) %>%
  rename_with(~ str_remove(., "^Fenton_"), .cols = everything()) %>%
  select(-Ontology)

df_H2O2 <- cmbd_lipids_H2O2 %>%
  mutate(Metabolite.curated = paste0(Metabolite.curated, " (H2O2)")) %>%
  rename_with(~ str_remove(., "^H2O2_"), .cols = everything()) %>%
  select(-Ontology)

df_stacked <- bind_rows(df_Fenton, df_H2O2)

# 2.3 Keep only oxidized lipids
cmbd_lipids <- df_stacked %>%
  filter(grepl("<", Metabolite.curated, fixed = TRUE))

df <- cmbd_lipids %>%
  ungroup() %>%
  column_to_rownames(var = colnames(.)[1])

# ============================================================
# 3. Reshape to long format and join metadata # -----
# ============================================================

# 3.1 Strip Fenton_/H2O2_ prefixes from Sample_info to match renamed column names
Sample_info_updated <- Sample_info %>%
  mutate(Sample_name_clean = str_remove(Sample_name, "^(Fenton_|H2O2_)"))

df_long <- df %>%
  rownames_to_column(var = "Lipid") %>%
  pivot_longer(-Lipid, names_to = "Sample", values_to = "Expression") %>%
  left_join(Sample_info_updated %>% select(Sample_name_clean, Time),
            by = c("Sample" = "Sample_name_clean")) %>%
  filter(!is.na(Time) & Time != "")

# 3.2 Add treatment column derived from lipid name suffix
df_long <- df_long %>%
  mutate(Treatment = ifelse(grepl("\\(Fenton\\)", Lipid), "Fenton", "H2O2"))

# ============================================================
# 4. Z-score normalization and K-means input matrix # -----
# ============================================================

# 4.1 Z-score normalization
df_zscore <- df_long %>%
  group_by(Lipid) %>%
  mutate(
    z_mean  = mean(Expression, na.rm = TRUE),
    z_sd    = sd(Expression,   na.rm = TRUE),
    Z_score = (Expression - z_mean) / z_sd
  ) %>%
  ungroup()

# 4.2 Average Z-scores per lipid per time point and pivot to wide format
df_scaled_input <- df_zscore %>%
  group_by(Lipid, Time) %>%
  summarise(Z_mean = mean(Z_score, na.rm = TRUE), .groups = "drop") %>%
  mutate(Time = factor(Time, levels = c("0 h", "1 min", "20 min", "4 h", "24 h", "72 h"))) %>%
  pivot_wider(names_from = Time, values_from = Z_mean) %>%
  column_to_rownames("Lipid")

# 4.3 Remove incomplete and zero-variance rows
df_scaled_input <- df_scaled_input[complete.cases(df_scaled_input), ]

zero_var_rows   <- apply(df_scaled_input, 1, function(x) var(x, na.rm = TRUE) == 0 | is.na(var(x, na.rm = TRUE)))
df_scaled_input <- df_scaled_input[!zero_var_rows, ]

# ============================================================
# 5. Cluster-wise bootstrap stability across candidate k # -----
# ============================================================

# 5.1 Run bootstrap resampling for each candidate k (Hennig, C. (2007). Cluster-wise
# assessment of cluster stability. Computational Statistics & Data Analysis, 52(1),
# 258-271. doi:10.1016/j.csda.2006.11.025).

cb_results <- lapply(K_SCAN, function(kk) {
  message("Bootstrapping k = ", kk, " (B = ", BOOT_B, ") ...")
  clusterboot(df_scaled_input,
              B             = BOOT_B,
              bootmethod    = "boot",
              clustermethod = kmeansCBI,
              krange        = kk,
              runs          = KM_NSTART,
              seed          = SEED,
              count         = FALSE)
})
names(cb_results) <- as.character(K_SCAN)

stability_scan <- lapply(K_SCAN, function(kk) {
  cb <- cb_results[[as.character(kk)]]
  data.frame(
    k            = kk,
    Cluster      = 1:kk,
    Mean_Jaccard = round(cb$bootmean, 3),
    N_dissolved  = cb$bootbrd,
    N_recovered  = cb$bootrecover
  )
}) %>%
  bind_rows() %>%
  mutate(
    Prop_dissolved = round(N_dissolved / BOOT_B, 3),
    Prop_recovered = round(N_recovered / BOOT_B, 3),
    Interpretation = case_when(
      Mean_Jaccard > 0.85 ~ "Highly stable",
      Mean_Jaccard > 0.75 ~ "Stable",
      Mean_Jaccard > 0.60 ~ "Pattern indicated, not reliable",
      TRUE                ~ "Not a valid cluster"
    )
  )

# 5.2 Summarise across k
scan_summary <- stability_scan %>%
  group_by(k) %>%
  summarise(
    n_stable     = sum(Mean_Jaccard > 0.75),
    n_total      = n(),
    prop_stable  = round(n_stable / n_total, 2),
    n_borderline = sum(Mean_Jaccard > 0.60 & Mean_Jaccard <= 0.75),
    n_invalid    = sum(Mean_Jaccard <= 0.60),
    min_jaccard  = min(Mean_Jaccard),
    mean_jaccard = round(mean(Mean_Jaccard), 3),
    sd_jaccard   = round(sd(Mean_Jaccard), 3),
    .groups = "drop"
  )

print(stability_scan)
print(as.data.frame(scan_summary))

write.csv(stability_scan, "./Output/cluster_stability_scan.csv", row.names = FALSE)
write.csv(scan_summary,   "./Output/cluster_stability_summary.csv", row.names = FALSE)

# 5.3 Stability profile plot across k
p_scan <- ggplot(stability_scan, aes(x = factor(Cluster), y = Mean_Jaccard)) +
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "#E05A51") +
  geom_hline(yintercept = 0.60, linetype = "dotted",  color = "grey40") +
  geom_col(aes(fill = Mean_Jaccard > 0.75), width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "#228B22", "FALSE" = "grey70"),
                    guide = "none") +
  facet_wrap(~ k, nrow = 2, scales = "free_x",
             labeller = labeller(k = function(x) paste0("k = ", x))) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Cluster", y = "Mean Jaccard similarity",
       title    = "Cluster-wise bootstrap stability",
       subtitle = paste0("B = ", BOOT_B,
                         "; dashed line = 0.75 (stable), dotted line = 0.60")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "grey95"))

print(p_scan)

ggsave("./Figure/cluster_stability_scan.svg", p_scan,
       width = 9, height = 5.5, bg = "transparent")

# 5.4 Mean stability across k
p_mean <- ggplot(scan_summary, aes(x = k, y = mean_jaccard)) +
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "#E05A51") +
  geom_line(color = "grey40") +
  geom_point(aes(color = k == K_CLUSTERS), size = 2.5) +
  scale_color_manual(values = c("TRUE" = "#E05A51", "FALSE" = "grey60"),
                     guide = "none") +
  scale_x_continuous(breaks = K_SCAN) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Number of clusters (k)", y = "Mean Jaccard similarity",
       title    = "Mean cluster stability across candidate k",
       subtitle = paste0("Red point = k used for the final figure (k = ", K_CLUSTERS, ")")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold", size = 11))

print(p_mean)

ggsave("./Figure/cluster_stability_mean.svg", p_mean,
       width = 6, height = 3.5, bg = "transparent")

# ============================================================
# 6. Final clustering at the selected k # -----
# ============================================================

k_use <- K_CLUSTERS

cb_final           <- cb_results[[as.character(k_use)]]
cluster_assignment <- cb_final$partition   # named vector, labels match stability_scan

cluster_stability <- stability_scan %>% filter(k == k_use)
print(cluster_stability)

message(sprintf("k = %d: mean (SD) Jaccard similarity = %.2f (%.2f) across all clusters",
                k_use,
                mean(cluster_stability$Mean_Jaccard),
                sd(cluster_stability$Mean_Jaccard)))

k <- k_use   # kept for the plotting code below

# ============================================================
# 7. Prepare scaled long format for plotting # -----
# ============================================================

df_scaled_long <- df_scaled_input %>%
  as.data.frame() %>%
  rownames_to_column(var = "Lipid") %>%
  pivot_longer(-Lipid, names_to = "Time", values_to = "Scaled_Expression") %>%
  mutate(Time = factor(Time, levels = c("0 h", "1 min", "20 min", "4 h", "24 h", "72 h")))

df_scaled_long <- df_scaled_long %>%
  mutate(Treatment = ifelse(grepl("\\(Fenton\\)", Lipid), "Fenton", "H2O2"))

df_clustered <- df_scaled_long %>%
  left_join(data.frame(Lipid   = names(cluster_assignment),
                       Cluster = as.integer(cluster_assignment)), by = "Lipid")

df_clustered$Cluster <- factor(df_clustered$Cluster, levels = 1:k)

n_lipids <- df_clustered %>%
  group_by(Cluster, Treatment) %>%
  summarise(n = n_distinct(Lipid), .groups = "drop")

max_n         <- max(n_lipids$n)
min_linewidth <- 0.5
max_linewidth <- 2.2

# ============================================================
# 8. Build one plot per cluster # -----
# ============================================================

plot_list <- lapply(1:k, function(i) {
  cluster_data <- df_clustered %>% filter(Cluster == i)
  
  n_fenton <- n_lipids %>% filter(Cluster == i, Treatment == "Fenton") %>% pull(n)
  n_h2o2   <- n_lipids %>% filter(Cluster == i, Treatment == "H2O2")   %>% pull(n)
  
  if (length(n_fenton) == 0) n_fenton <- 0
  if (length(n_h2o2)   == 0) n_h2o2   <- 0
  
  cluster_mean <- lapply(c("Fenton", "H2O2"), function(trt) {
    lipids_in_group <- cluster_data %>%
      filter(Treatment == trt) %>%
      pull(Lipid) %>%
      unique()
    if (length(lipids_in_group) == 0) return(NULL)
    df_scaled_long %>%
      filter(Lipid %in% lipids_in_group) %>%
      group_by(Time) %>%
      summarise(
        Cluster_Mean = mean(Scaled_Expression, na.rm = TRUE),
        Treatment    = trt,
        .groups = "drop"
      )
  }) %>% bind_rows()
  
  cluster_mean <- cluster_mean %>%
    mutate(
      n_lipids_treatment = ifelse(Treatment == "Fenton", n_fenton, n_h2o2),
      linewidth_scaled   = min_linewidth + (n_lipids_treatment / max_n) *
        (max_linewidth - min_linewidth)
    )
  
  y_range  <- range(cluster_data$Scaled_Expression, na.rm = TRUE)
  y_limits <- c(y_range[1] - 0.5, y_range[2] + 0.5)
  
  ggplot(cluster_data, aes(x = Time, y = Scaled_Expression, group = Lipid)) +
    geom_line(alpha = 0.3, linewidth = 0.5, color = "grey60") +
    geom_line(data = cluster_mean,
              aes(x = Time, y = Cluster_Mean, group = Treatment,
                  color = Treatment, linewidth = linewidth_scaled),
              alpha = 1) +
    scale_linewidth_identity() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    scale_color_manual(values = c("Fenton" = "#E05A51", "H2O2" = "#377EB8")) +
    coord_cartesian(ylim = y_limits) +
    labs(title    = bquote(bold(.(paste("Cluster", i)))),
         subtitle = bquote("(Fenton:"~.(n_fenton)*","~H[2]*O[2]*"-only:"~.(n_h2o2)*")"),
         y = "Z score") +
    theme_bw() +
    theme(
      legend.position  = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(t = 9, r = 10, b = 3, l = 3, unit = "pt"),
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 11,
                                      margin = margin(t = 2, b = 2)),
      plot.subtitle    = element_text(hjust = 0.5, size = 10,
                                      margin = margin(t = 3, b = 3)),
      axis.title       = element_text(size = 8),
      axis.title.x     = element_blank(),
      axis.title.y     = element_text(size = 10),
      axis.text        = element_text(size = 8.5),
      axis.text.x      = element_text(angle = 45, hjust = 1)
    )
})

# ============================================================
# 9. Shared legend # -----
# ============================================================

dummy_data <- data.frame(
  x         = c(1, 2),
  xend      = c(2, 3),
  y         = c(1, 2),
  yend      = c(1, 2),
  Treatment = factor(c("Fenton", "H2O2"), levels = c("Fenton", "H2O2"))
)

dummy_plot <- ggplot(dummy_data, aes(x = x, xend = xend,
                                     y = y, yend = yend,
                                     color = Treatment)) +
  geom_segment(linewidth = 2.5) +
  scale_color_manual(
    name   = "Oxidized lipid references",
    values = c("Fenton" = "#E05A51", "H2O2" = "#377EB8"),
    labels = c("Fenton" = "Fenton", "H2O2" = expression(H[2]*O[2]*"-only"))
  ) +
  theme_bw() +
  theme(
    legend.position   = "bottom",
    legend.title      = element_text(size = 11, face = "bold"),
    legend.text       = element_text(size = 10),
    legend.key.width  = unit(0.6, "cm"),
    legend.key.height = unit(0.1, "cm")
  )

legend <- get_legend(dummy_plot)

# ============================================================
# 10. Assemble, export, and save results # -----
# ============================================================

# 10.1 Column spacing control
col_gap  <- 0.1   # spacer width relative to one plot column (0 = no gap)

# 10.2 Assemble all clusters in a single row, with spacers between columns
spacer <- nullGrob()
n_rows <- 1
n_cols <- k

grob_list <- list()
for (i in seq_len(n_cols)) {
  grob_list <- c(grob_list, list(plot_list[[i]]))
  if (i < n_cols) grob_list <- c(grob_list, list(spacer))
}

widths_vec <- rep(1, 2 * n_cols - 1)
widths_vec[seq(2, length(widths_vec), by = 2)] <- col_gap

grid_body <- arrangeGrob(
  grobs  = grob_list,
  ncol   = 2 * n_cols - 1,
  nrow   = n_rows,
  widths = unit(widths_vec, "null")
)

final_plot <- arrangeGrob(grid_body, legend, ncol = 1,
                          heights = c(n_rows * (10 / 3), 0.5))

# 10.3 Preview
grid.newpage()
grid.draw(final_plot)

# 10.4 Export as SVG with width scaled to the number of columns
svglite(paste0("./Figure/cluster_trends_Fenton_H2O2_k", k, ".svg"),
        width  = (n_cols * 200) * 1 / 96,
        height = (n_rows * 200 + 30) * 1 / 96,
        bg     = "transparent")
grid.draw(final_plot)
dev.off()

# 10.5 Save clustering results
Fenton_H2O2_clusters <- df_clustered %>%
  select(Lipid, Cluster) %>%
  distinct()

save(Fenton_H2O2_clusters, file = "./Output/Fenton_H2O2_clusters.Rdata")
save(stability_scan, scan_summary, cluster_stability,
     file = "./Output/Fenton_H2O2_cluster_stability.Rdata")