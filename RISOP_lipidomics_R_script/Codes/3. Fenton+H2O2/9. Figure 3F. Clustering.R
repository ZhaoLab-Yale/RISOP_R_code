rm(list = ls())
library(tidyverse)
library(gridExtra)
library(grid)
library(ggpubr)
library(svglite)
setwd(rstudioapi::getActiveProject())

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

# 3.2 Add replicate index
df_long <- df_long %>%
  mutate(
    Replicate = case_when(
      grepl("_1$|_Rep1$|_r1$", Sample, ignore.case = TRUE) ~ "1",
      grepl("_2$|_Rep2$|_r2$", Sample, ignore.case = TRUE) ~ "2",
      grepl("_3$|_Rep3$|_r3$", Sample, ignore.case = TRUE) ~ "3",
      TRUE ~ as.character((row_number() - 1) %% 3 + 1)
    )
  )

# 3.3 Add treatment column derived from lipid name suffix
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

# 4.3 Remove zero-variance rows
zero_var_rows   <- apply(df_scaled_input, 1, function(x) var(x, na.rm = TRUE) == 0 | is.na(var(x, na.rm = TRUE)))
df_scaled_input <- df_scaled_input[!zero_var_rows, ]

if (any(zero_var_rows))
  message("Zero-variance rows removed before K-means: ", sum(zero_var_rows))

# ============================================================
# 5. K-means clustering # -----
# ============================================================

k <- min(9, nrow(df_scaled_input) - 1)
message("Running K-means with k = ", k)

set.seed(123)
kmeans_result <- kmeans(df_scaled_input, centers = k, nstart = 25)

# ============================================================
# 6. Prepare scaled long format for plotting # -----
# ============================================================

df_scaled_long <- df_scaled_input %>%
  as.data.frame() %>%
  rownames_to_column(var = "Lipid") %>%
  pivot_longer(-Lipid, names_to = "Time", values_to = "Scaled_Expression") %>%
  mutate(Time = factor(Time, levels = c("0 h", "1 min", "20 min", "4 h", "24 h", "72 h")))

df_scaled_long <- df_scaled_long %>%
  mutate(Treatment = ifelse(grepl("\\(Fenton\\)", Lipid), "Fenton", "H2O2"))

df_clustered <- df_scaled_long %>%
  left_join(data.frame(Lipid   = names(kmeans_result$cluster),
                       Cluster = kmeans_result$cluster), by = "Lipid")

df_clustered$Cluster <- factor(df_clustered$Cluster, levels = 1:k)

n_lipids <- df_clustered %>%
  group_by(Cluster, Treatment) %>%
  summarise(n = n_distinct(Lipid), .groups = "drop")

max_n         <- max(n_lipids$n)
min_linewidth <- 0.8
max_linewidth <- 3

# ============================================================
# 7. Build one plot per cluster # -----
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
         x = "Time",
         y = "Z score") +
    theme_bw() +
    theme(
      legend.position  = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(t = 9, r = 10, b = 3, l = 8, unit = "pt"),
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 11,
                                      margin = margin(t = 2, b = 2)),
      plot.subtitle    = element_text(hjust = 0.5, size = 10,
                                      margin = margin(t = 3, b = 3)),
      axis.title       = element_text(size = 8),
      axis.title.x     = element_text(size = 10),
      axis.title.y     = element_text(size = 10),
      axis.text        = element_text(size = 8.5),
      axis.text.x      = element_text(angle = 45, hjust = 1)
    )
})

# ============================================================
# 8. Shared legend # -----
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
# 9. Assemble, export, and save results # -----
# ============================================================

# 9.1 Assemble final plot
final_plot <- arrangeGrob(
  arrangeGrob(grobs = plot_list, ncol = 3, nrow = 3),
  legend,
  ncol    = 1,
  heights = c(10, 0.5)
)

# 9.2 Preview
grid.arrange(
  arrangeGrob(grobs = plot_list, ncol = 3, nrow = 3),
  legend,
  ncol    = 1,
  heights = c(10, 0.5)
)

# 9.3 Export as SVG
svglite("./Figure/cluster_trends_Fenton_H2O2.svg",
        width  = 650 / 96,
        height = 700 / 96,
        bg     = "transparent")
grid.draw(final_plot)
dev.off()

# 9.4 Save clustering results
Fenton_H2O2_clusters <- df_clustered %>%
  select(Lipid, Cluster) %>%
  distinct()

save(Fenton_H2O2_clusters, file = "./Output/Fenton_H2O2_clusters.Rdata")