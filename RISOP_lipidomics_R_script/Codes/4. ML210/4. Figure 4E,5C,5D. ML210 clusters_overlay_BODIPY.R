rm(list = ls())
library(tidyverse)
library(gridExtra)
library(grid)
library(egg)
library(ggtext)
setwd(rstudioapi::getActiveProject())

# ============================================================
# 1. Load data # -----
# ============================================================

load(file = "./Output/cmbd_lipids.Rdata")
Sample_info <- read.csv('./Input/Sample_info.csv', stringsAsFactors = F, skip = 0, header = T)
BODIPY      <- read.csv('./Input/BODIPY.csv',      stringsAsFactors = F, skip = 0, header = T)

select <- dplyr::select
filter <- dplyr::filter

# ============================================================
# 2. Shared time-point levels and cluster colors # -----
# ============================================================

time_levels_numeric  <- c("0", "20", "45", "80", "120", "270")
time_levels_original <- c("0 min", "20 min", "45 min", "80 min", "120 min", "270 min")

cluster_colors <- c("1" = "#A1863B",
                    "2" = "#BD8349",
                    "3" = "#8B5742")

# ============================================================
# 3. Prepare OxPL expression matrix # -----
# ============================================================

# 3.1 Filter to oxidized lipids
cmbd_lipids <- cmbd_lipids %>%
  dplyr::filter(grepl("<", Metabolite.curated, fixed = TRUE))

message("Oxidized lipids: ", nrow(cmbd_lipids))

# 3.2 Pivot to long format and join time metadata
df <- cmbd_lipids %>%
  ungroup() %>%
  dplyr::select(-Ontology) %>%
  column_to_rownames(var = colnames(.)[1])

df_long <- df %>%
  rownames_to_column(var = "Lipid") %>%
  pivot_longer(-Lipid, names_to = "Sample", values_to = "Expression") %>%
  left_join(Sample_info %>% dplyr::select(Sample_name, Time),
            by = c("Sample" = "Sample_name")) %>%
  dplyr::filter(!is.na(Time) & Time != "")

# 3.3 Add replicate index
df_long <- df_long %>%
  mutate(
    Replicate = case_when(
      grepl("_1$|_Rep1$|_r1$", Sample, ignore.case = TRUE) ~ "1",
      grepl("_2$|_Rep2$|_r2$", Sample, ignore.case = TRUE) ~ "2",
      grepl("_3$|_Rep3$|_r3$", Sample, ignore.case = TRUE) ~ "3",
      TRUE ~ as.character((row_number() - 1) %% 3 + 1)
    )
  )

# ============================================================
# 4. Z-score normalization and K-means clustering # -----
# ============================================================

# 4.1 Z-score normalize per lipid
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
  mutate(Time = factor(Time, levels = time_levels_original,
                       labels = time_levels_numeric)) %>%
  pivot_wider(names_from = Time, values_from = Z_mean) %>%
  column_to_rownames("Lipid")

# 4.3 Remove zero-variance rows
zero_var_rows   <- apply(df_scaled_input, 1,
                         function(x) var(x, na.rm = TRUE) == 0 | is.na(var(x, na.rm = TRUE)))
df_scaled_input <- df_scaled_input[!zero_var_rows, ]

if (any(zero_var_rows))
  message("Zero-variance rows removed before K-means: ", sum(zero_var_rows))

# 4.4 Run K-means (k = 3)
set.seed(12345678)
k             <- 3
kmeans_result <- kmeans(df_scaled_input, centers = k, nstart = 25)

# ============================================================
# 5. Prepare scaled long format and cluster assignments # -----
# ============================================================

df_scaled_long <- df_scaled_input %>%
  as.data.frame() %>%
  rownames_to_column(var = "Lipid") %>%
  pivot_longer(-Lipid, names_to = "Time", values_to = "Scaled_Expression") %>%
  mutate(Time = factor(Time, levels = time_levels_numeric))

df_clustered <- df_scaled_long %>%
  left_join(data.frame(Lipid   = names(kmeans_result$cluster),
                       Cluster = kmeans_result$cluster), by = "Lipid")

df_clustered$Cluster <- factor(df_clustered$Cluster, levels = 1:k)

n_lipids <- df_clustered %>%
  group_by(Cluster) %>%
  summarise(n = n_distinct(Lipid))

max_n         <- max(n_lipids$n)
min_linewidth <- 0.5
max_linewidth <- 1.5

# ============================================================
# 6. Compute cluster mean curves # -----
# ============================================================

all_cluster_means <- lapply(1:k, function(i) {
  lipids_in_cluster <- df_clustered %>%
    dplyr::filter(Cluster == i) %>%
    pull(Lipid) %>%
    unique()
  
  df_scaled_long %>%
    dplyr::filter(Lipid %in% lipids_in_cluster) %>%
    group_by(Time) %>%
    summarise(
      Cluster_Mean = mean(Scaled_Expression, na.rm = TRUE),
      SEM          = sd(Scaled_Expression,   na.rm = TRUE) / sqrt(n_distinct(Lipid)),
      Cluster      = as.factor(i),
      .groups = "drop"
    )
}) %>% bind_rows()

all_cluster_means <- all_cluster_means %>%
  left_join(n_lipids, by = "Cluster") %>%
  mutate(linewidth_scaled = min_linewidth + (n / max_n) * (max_linewidth - min_linewidth))

# ============================================================
# 7. Cluster line plots # -----
# ============================================================

plot_list <- lapply(1:k, function(i) {
  cluster_data     <- df_clustered %>% dplyr::filter(Cluster == i)
  cluster_mean     <- all_cluster_means %>% dplyr::filter(Cluster == i)
  n                <- n_lipids %>% dplyr::filter(Cluster == i) %>% pull(n)
  linewidth_scaled <- unique(cluster_mean$linewidth_scaled)
  line_color       <- cluster_colors[as.character(i)]
  
  y_range  <- range(cluster_data$Scaled_Expression, na.rm = TRUE)
  y_limits <- c(y_range[1] - 0.5, y_range[2] + 0.5)
  
  ggplot(cluster_data, aes(x = Time, y = Scaled_Expression, group = Lipid)) +
    geom_line(alpha = 0.4, linewidth = 0.5, color = "gray60") +
    geom_line(data = cluster_mean,
              aes(x = Time, y = Cluster_Mean, group = 1),
              color = line_color, linewidth = linewidth_scaled) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    coord_cartesian(ylim = y_limits) +
    labs(title = paste0("Cluster ", i, " (n = ", n, ")"),
         x = "Time (min)", y = "Z score") +
    theme_bw() +
    theme(text             = element_text(family = "Helvetica"),
          legend.position  = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          plot.title       = element_text(hjust = 0.5, face = "bold", size = 11),
          axis.title       = element_text(size = 10, face = "bold"),
          axis.title.x     = element_text(size = 10, face = "bold", vjust = -0.5),
          axis.text        = element_text(size = 8),
          axis.text.y      = element_text(size = 8.5),
          axis.text.x      = element_text(angle = 0, hjust = 0.5, size = 8.5))
})

grid.arrange(grobs = plot_list, ncol = 1, nrow = 3)

# ============================================================
# 8. Prepare BODIPY C11 data # -----
# ============================================================

# 8.1 Pivot to long format and scale
BODIPY_long <- BODIPY %>%
  pivot_longer(cols = -Replicate,
               names_to  = "Time_raw",
               values_to = "Value") %>%
  dplyr::filter(!is.na(Value)) %>%
  mutate(Time_raw2 = gsub("X",    "",   Time_raw),
         Time_raw2 = gsub("\\.", " ", Time_raw2),
         Time      = factor(Time_raw2,
                            levels = time_levels_original,
                            labels = time_levels_numeric),
         Replicate = as.character(Replicate))

mean_bodipy <- mean(BODIPY_long$Value, na.rm = TRUE)
sd_bodipy   <- sd(BODIPY_long$Value,   na.rm = TRUE)

BODIPY_scaled <- BODIPY_long %>%
  mutate(Scaled_Value = (Value - mean_bodipy) / sd_bodipy)

BODIPY_stats <- BODIPY_scaled %>%
  group_by(Time) %>%
  summarise(
    Mean_Scaled = mean(Scaled_Value, na.rm = TRUE),
    SEM         = sd(Scaled_Value,   na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# 8.2 Normality check (Shapiro-Wilk)
normality_results <- BODIPY_long %>%
  group_by(Time) %>%
  summarise(
    n         = n(),
    statistic = shapiro.test(Value)$statistic,
    p_value   = shapiro.test(Value)$p.value,
    normal    = ifelse(shapiro.test(Value)$p.value > 0.05, "Yes", "No"),
    .groups   = "drop"
  )

message("=== Shapiro-Wilk Normality Test (BODIPY per time point) ===")
print(normality_results)

normality_qq_plot <- ggplot(BODIPY_long, aes(sample = Value)) +
  stat_qq(color = "#228B22", size = 2) +
  stat_qq_line(color = "black", linetype = "dashed") +
  facet_wrap(~ Time, ncol = 3,
             labeller = labeller(Time = function(x) paste0(x, " min"))) +
  labs(title = "Q-Q Plots: BODIPY C11 by Time Point",
       x     = "Theoretical Quantiles",
       y     = "Sample Quantiles") +
  theme_bw() +
  theme(text             = element_text(family = "Helvetica"),
        plot.title       = element_text(hjust = 0.5, face = "bold", size = 12),
        strip.background = element_rect(fill = "grey90"),
        strip.text       = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

print(normality_qq_plot)

# ============================================================
# 9. Compute cluster total intensity (sum then Z-score) # -----
# ============================================================

# 9.1 Map cluster labels onto raw expression data
lipid_cluster_map <- data.frame(
  Lipid   = names(kmeans_result$cluster),
  Cluster = as.factor(kmeans_result$cluster)
)

df_long_clustered <- df_long %>%
  left_join(lipid_cluster_map, by = "Lipid") %>%
  dplyr::filter(!is.na(Cluster))

# 9.2 Sum all lipid intensities per cluster per sample
df_cluster_total <- df_long_clustered %>%
  group_by(Cluster, Sample, Time) %>%
  summarise(Total_Intensity = sum(Expression, na.rm = TRUE), .groups = "drop")

# 9.3 Z-score within each cluster
df_cluster_total <- df_cluster_total %>%
  group_by(Cluster) %>%
  mutate(
    mean_total       = mean(Total_Intensity, na.rm = TRUE),
    sd_total         = sd(Total_Intensity,   na.rm = TRUE),
    Scaled_Intensity = (Total_Intensity - mean_total) / sd_total
  ) %>%
  ungroup()

# 9.4 Summarize mean and SEM per cluster per time point
cluster_total_summary <- df_cluster_total %>%
  mutate(Time = factor(Time, levels = time_levels_original,
                       labels = time_levels_numeric)) %>%
  group_by(Cluster, Time) %>%
  summarise(
    Mean_Total_Scaled = mean(Scaled_Intensity, na.rm = TRUE),
    SEM_Total_Scaled  = sd(Scaled_Intensity,   na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# ============================================================
# 10. Pearson correlations: cluster total Z-score vs BODIPY # -----
# ============================================================

correlation_data <- BODIPY_stats %>%
  dplyr::select(Time, BODIPY_Mean = Mean_Scaled) %>%
  left_join(
    cluster_total_summary %>%
      dplyr::select(Time, Cluster, Mean_Total_Scaled) %>%
      pivot_wider(names_from = Cluster, values_from = Mean_Total_Scaled,
                  names_prefix = "Cluster_"),
    by = "Time"
  )

correlations <- data.frame(
  Comparison = character(),
  Cluster    = character(),
  Pearson_r  = numeric(),
  Pearson_p  = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:k) {
  cluster_col  <- paste0("Cluster_", i)
  pearson_test <- cor.test(correlation_data$BODIPY_Mean,
                           correlation_data[[cluster_col]],
                           method = "pearson")
  correlations <- rbind(correlations, data.frame(
    Comparison = paste0("BODIPY vs Cluster ", i),
    Cluster    = as.character(i),
    Pearson_r  = pearson_test$estimate,
    Pearson_p  = pearson_test$p.value
  ))
}

correlations$Pearson_p_adj <- p.adjust(correlations$Pearson_p, method = "BH")
correlations <- correlations %>%
  mutate(
    Pearson_sig = case_when(
      Pearson_p_adj < 0.001 ~ "***",
      Pearson_p_adj < 0.01  ~ "**",
      Pearson_p_adj < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

print(correlations)

# ============================================================
# 11. Build overlay plot (cluster total Z-score vs BODIPY) # -----
# ============================================================

# 11.1 Legend labels with correlation statistics
legend_labels <- c(
  sprintf("Total OxPLs in cluster 1 (*R* = %.3f, *P* = %.3f)",
          correlations$Pearson_r[correlations$Cluster == "1"],
          correlations$Pearson_p_adj[correlations$Cluster == "1"]),
  sprintf("Total OxPLs in cluster 2 (*R* = %.3f, *P* = %.3f)",
          correlations$Pearson_r[correlations$Cluster == "2"],
          correlations$Pearson_p_adj[correlations$Cluster == "2"]),
  sprintf("Total OxPLs in cluster 3 (*R* = %.3f, *P* = %.3f)",
          correlations$Pearson_r[correlations$Cluster == "3"],
          correlations$Pearson_p_adj[correlations$Cluster == "3"]),
  "Peroxidation\u2212positive cells by BODIPY C11"
)

# 11.2 Bracket segments showing significance at last time point
last_time <- "270"

endpoint_values <- bind_rows(
  cluster_total_summary %>%
    dplyr::filter(Time == last_time) %>%
    dplyr::select(Cluster, y_end = Mean_Total_Scaled) %>%
    mutate(Cluster = as.character(Cluster)),
  BODIPY_stats %>%
    dplyr::filter(Time == last_time) %>%
    dplyr::select(y_end = Mean_Scaled) %>%
    mutate(Cluster = "BODIPY")
)

bracket_data <- data.frame(Cluster = as.character(1:k)) %>%
  left_join(
    endpoint_values %>% dplyr::filter(Cluster != "BODIPY") %>% rename(y_cluster = y_end),
    by = "Cluster"
  ) %>%
  mutate(
    y_bodipy = endpoint_values %>% dplyr::filter(Cluster == "BODIPY") %>% pull(y_end),
    y_min    = pmin(y_cluster, y_bodipy),
    y_max    = pmax(y_cluster, y_bodipy),
    y_mid    = (y_min + y_max) / 2
  ) %>%
  left_join(correlations %>% dplyr::select(Cluster, Pearson_sig), by = "Cluster") %>%
  rename(Significance = Pearson_sig) %>%
  mutate(
    text_size = case_when(
      Significance == "ns" ~ 3.5,
      TRUE ~ 6
    )
  )

x_bracket_base    <- 6.45
bracket_widths    <- c(0, 0.3, 1.0)
bracket_linewidth <- 0.3

bracket_segments <- data.frame()
for (i in 1:nrow(bracket_data)) {
  row      <- bracket_data[i, ]
  x_offset <- bracket_widths[i]
  bracket_segments <- rbind(
    bracket_segments,
    data.frame(Cluster = row$Cluster,
               x = x_bracket_base + x_offset, xend = x_bracket_base + x_offset,
               y = row$y_min, yend = row$y_max, segment_type = "vertical"),
    data.frame(Cluster = row$Cluster,
               x = x_bracket_base + x_offset - 0.2, xend = x_bracket_base + x_offset,
               y = row$y_max, yend = row$y_max, segment_type = "top"),
    data.frame(Cluster = row$Cluster,
               x = x_bracket_base + x_offset - 0.2, xend = x_bracket_base + x_offset,
               y = row$y_min, yend = row$y_min, segment_type = "bottom")
  )
}

bracket_data <- bracket_data %>%
  mutate(x_label = x_bracket_base + bracket_widths + 0.1)

# 11.3 Assemble overlay plot
overlay_plot_final <- ggplot() +
  geom_ribbon(data = cluster_total_summary,
              aes(x = Time, y = Mean_Total_Scaled,
                  ymin = Mean_Total_Scaled - SEM_Total_Scaled,
                  ymax = Mean_Total_Scaled + SEM_Total_Scaled,
                  group = Cluster, fill = Cluster),
              alpha = 0.2, color = NA) +
  geom_line(data = cluster_total_summary,
            aes(x = Time, y = Mean_Total_Scaled,
                group = Cluster, color = Cluster),
            linewidth = 1.0) +
  geom_point(data = cluster_total_summary,
             aes(x = Time, y = Mean_Total_Scaled,
                 color = Cluster, shape = Cluster),
             size = 3) +
  geom_ribbon(data = BODIPY_stats,
              aes(x = Time, y = Mean_Scaled,
                  ymin = Mean_Scaled - SEM, ymax = Mean_Scaled + SEM,
                  group = 1, fill = "BODIPY"),
              alpha = 0.2, color = NA) +
  geom_line(data = BODIPY_stats,
            aes(x = Time, y = Mean_Scaled, group = 1, color = "BODIPY"),
            linewidth = 1.0) +
  geom_point(data = BODIPY_stats,
             aes(x = Time, y = Mean_Scaled, color = "BODIPY", shape = "BODIPY"),
             size = 3.5) +
  geom_segment(data = bracket_segments,
               aes(x = x, xend = xend, y = y, yend = yend),
               color = "black", linewidth = bracket_linewidth) +
  geom_text(data = bracket_data,
            aes(x = x_label, y = y_mid, label = Significance, size = text_size),
            color = "black", hjust = 0, vjust = 0.65, fontface = "bold",
            family = "Helvetica") +
  scale_size_identity() +
  scale_color_manual(values = c(cluster_colors, "BODIPY" = "#228B22"),
                     labels = legend_labels, name = "") +
  scale_fill_manual(values  = c(cluster_colors, "BODIPY" = "#228B22"),
                    labels  = legend_labels, name = "") +
  scale_shape_manual(values = c("1" = 17, "2" = 17, "3" = 17, "BODIPY" = 16),
                     labels = legend_labels, name = "") +
  guides(shape = guide_legend(override.aes = list(size = c(3, 3, 3, 3)),
                              keyheight = unit(1.3, "lines")),
         color = guide_legend(override.aes = list(size = c(3, 3, 3, 3)),
                              keyheight = unit(1.3, "lines"))) +
  labs(x = "ML210 treatment time (min)", y = "Z score") +
  coord_cartesian(clip = "off", xlim = c(0.8, 7.5), ylim = c(-1.9, 2.4)) +
  theme_bw() +
  theme(text                 = element_text(family = "Helvetica"),
        panel.grid.minor     = element_blank(),
        panel.grid.major     = element_blank(),
        panel.border         = element_rect(color = "black", linewidth = 0.3),
        plot.title           = element_markdown(hjust = 0.5, size = 12,
                                                family = "Helvetica"),
        axis.title           = element_text(size = 13.05),
        axis.title.x         = element_text(size = 13.05, vjust = -1.5),
        axis.text            = element_text(size = 11.5),
        axis.text.x          = element_text(angle = 0, hjust = 0.5),
        legend.position      = "right",
        legend.justification = c(0, 0),
        legend.box.margin    = margin(0, 0, 0, 0),
        legend.title         = element_text(face = "bold"),
        legend.text          = element_markdown(size = 9, family = "Helvetica"),
        aspect.ratio         = 0.75,
        plot.margin          = unit(c(0, 1, 0, 3), "mm"))

print(overlay_plot_final)

# ============================================================
# 12. BODIPY C11 scatter plot with significance # -----
# ============================================================

# 12.1 Compute summary statistics
summary_stats <- BODIPY_long %>%
  group_by(Time) %>%
  summarise(
    median = median(Value, na.rm = TRUE),
    mean   = mean(Value,   na.rm = TRUE),
    sd     = sd(Value,     na.rm = TRUE),
    n      = n(),
    se     = sd / sqrt(n)
  )

# 12.2 t-tests vs baseline with BH correction
baseline_data    <- BODIPY_long %>% dplyr::filter(Time == "0") %>% pull(Value)
time_levels_test <- c("20", "45", "80", "120", "270")

p_values   <- sapply(time_levels_test, function(t) {
  test_data <- BODIPY_long %>% dplyr::filter(Time == t) %>% pull(Value)
  t.test(test_data, baseline_data)$p.value
})
p_adjusted <- p.adjust(p_values, method = "BH")

sig_labels <- case_when(
  p_adjusted < 0.001 ~ "***",
  p_adjusted < 0.01  ~ "**",
  p_adjusted < 0.05  ~ "*",
  TRUE ~ "ns"
)

sig_data <- data.frame(
  Time         = factor(time_levels_test, levels = time_levels_numeric),
  p_value      = p_values,
  p_adjusted   = p_adjusted,
  significance = sig_labels
) %>%
  left_join(summary_stats, by = "Time") %>%
  mutate(y_position = mean + se + 5)

# 12.3 Build scatter plot
BODIPY_plot <- ggplot() +
  geom_point(data = BODIPY_long, aes(x = as.factor(Time), y = Value),
             position = position_jitter(width = 0.3, seed = 123),
             alpha = 0.3, size = 4, color = "#228B22") +
  geom_errorbar(data = summary_stats,
                aes(x = as.factor(Time), ymin = mean - se, ymax = mean + se),
                width = 0.3, linewidth = 0.8, color = "darkgreen") +
  geom_errorbar(data = summary_stats,
                aes(x = as.factor(Time), ymin = median, ymax = median),
                width = 0.5, linewidth = 0.3, color = "black") +
  geom_text(data = sig_data,
            aes(x = Time, y = y_position, label = significance),
            size = 7, vjust = 0.5, family = "Helvetica") +
  labs(x = "ML210 treatment time (min)", y = "Peroxidation-positive cells (%)") +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, by = 25)) +
  theme_bw() +
  theme(text             = element_text(family = "Helvetica"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title       = element_markdown(hjust = 0.5, size = 12,
                                            family = "Helvetica"),
        axis.title       = element_text(size = 12.5),
        axis.title.x     = element_text(size = 12.5, vjust = -2),
        axis.title.y     = element_text(margin = margin(r = 7)),
        axis.text        = element_text(size = 11.5),
        axis.text.x      = element_text(angle = 0, hjust = 0.5, size = 11.5),
        legend.position  = "right",
        legend.text      = element_text(size = 10),
        legend.title     = element_text(face = "bold"),
        aspect.ratio     = 0.75)

print(BODIPY_plot)

# ============================================================
# 13. Export figures and save clustering results # -----
# ============================================================

# 13.1 Export figures
pdf("./Figure/BODIPY_normality_QQ.pdf", width = 600 / 96, height = 400 / 96)
print(normality_qq_plot)
dev.off()

pdf("./Figure/clusters_plot.pdf", width = 600 / 96, height = 168 / 96)
grid.arrange(grobs = plot_list, ncol = 3, nrow = 1)
dev.off()

pdf("./Figure/overlay_plot.pdf", width = 645*1.07 / 96, height = 388*1.07 / 96)
print(overlay_plot_final)
dev.off()

pdf("./Figure/BODIPY_plot.pdf", width = 390 / 96, height = 330 / 96)
print(BODIPY_plot)
dev.off()

# 13.2 Save clustering results
ML210_clusters <- df_clustered %>%
  dplyr::select(Lipid, Cluster) %>%
  distinct()

save(ML210_clusters, file = "./Output/ML210_clusters.Rdata")
