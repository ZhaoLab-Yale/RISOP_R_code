rm(list = ls())
library(tidyverse)
library(openxlsx)
library(ggplot2)
library(ggtext)
setwd(rstudioapi::getActiveProject())

# ============================================================
# 1. Load data # -----
# ============================================================

load("./Output/cmbd_lipids.Rdata")
BODIPY         <- read.csv('./Input/BODIPY.csv',         stringsAsFactors = F, skip = 0, header = T)
species_for_QT <- read.csv('./Input/species_for_QT.csv', stringsAsFactors = F, skip = 0, header = T)
Sample_info    <- read.csv('./Input/Sample_info.csv',    stringsAsFactors = F, skip = 0, header = T)

# ============================================================
# 2. Shared time-point levels # -----
# ============================================================

time_order        <- c("0 min", "20 min", "45 min", "80 min", "120 min", "270 min")
time_labels_clean <- c("0", "20", "45", "80", "120", "270")
time_points       <- setdiff(time_order, "0 min")

# ============================================================
# 3. Prepare OxPL expression data # -----
# ============================================================

# 3.1 Filter to oxidized lipids
df_curated <- cmbd_lipids %>%
  filter(grepl("<", Metabolite.curated, fixed = TRUE))

# 3.2 Pivot to long format and join time metadata
data_long <- df_curated %>%
  pivot_longer(cols      = -c(Metabolite.curated, Ontology),
               names_to  = "Sample_name",
               values_to = "Intensity") %>%
  left_join(Sample_info %>% select(Sample_name, Time), by = "Sample_name") %>%
  filter(!is.na(Ontology) & !is.na(Time)) %>%
  mutate(Time = factor(Time, levels = time_order))

# ============================================================
# 4. Total OxPLs: sum and Z-score # -----
# ============================================================

# 4.1 Sum intensity per sample, then Z-score globally
data_total <- data_long %>%
  group_by(Sample_name, Time) %>%
  summarise(Total_Intensity = sum(Intensity, na.rm = TRUE), .groups = "drop")

mean_total_raw <- mean(data_total$Total_Intensity, na.rm = TRUE)
sd_total_raw   <- sd(data_total$Total_Intensity,   na.rm = TRUE)

data_total <- data_total %>%
  mutate(Scaled_Intensity = (Total_Intensity - mean_total_raw) / sd_total_raw)

# 4.2 Summarize mean and SEM per time point
data_total_summary <- data_total %>%
  group_by(Time) %>%
  summarise(
    Mean_Total_Scaled = mean(Scaled_Intensity, na.rm = TRUE),
    SEM_Total_Scaled  = sd(Scaled_Intensity,   na.rm = TRUE) / sqrt(n()),
    N                 = n(),
    .groups = "drop"
  )

# ============================================================
# 5. Normality tests and Q-Q plots # -----
# ============================================================

# 5.1 Prepare BODIPY long format for normality testing
BODIPY_long_for_norm <- BODIPY %>%
  pivot_longer(cols     = -Replicate,
               names_to  = "Time_raw",
               values_to = "Value") %>%
  filter(!is.na(Value)) %>%
  mutate(Time = gsub("X",    "", Time_raw),
         Time = gsub("\\.", " ", Time),
         Time = factor(Time, levels = time_order))

# 5.2 Shapiro-Wilk: Total OxPLs
normality_oxpl <- data_total %>%
  group_by(Time) %>%
  summarise(
    n         = n(),
    statistic = shapiro.test(Total_Intensity)$statistic,
    p_value   = shapiro.test(Total_Intensity)$p.value,
    normal    = ifelse(shapiro.test(Total_Intensity)$p.value > 0.05, "Yes", "No"),
    .groups   = "drop"
  )

message("=== Shapiro-Wilk Normality Test: Total OxPLs per time point ===")
print(normality_oxpl)

# 5.3 Shapiro-Wilk: BODIPY C11
normality_bodipy <- BODIPY_long_for_norm %>%
  group_by(Time) %>%
  summarise(
    n         = n(),
    statistic = shapiro.test(Value)$statistic,
    p_value   = shapiro.test(Value)$p.value,
    normal    = ifelse(shapiro.test(Value)$p.value > 0.05, "Yes", "No"),
    .groups   = "drop"
  )

message("=== Shapiro-Wilk Normality Test: BODIPY C11 per time point ===")
print(normality_bodipy)

# 5.4 Q-Q plot: Total OxPLs
qq_oxpl <- ggplot(data_total, aes(sample = Total_Intensity)) +
  stat_qq(color = "#7E1FF2", size = 2) +
  stat_qq_line(color = "black", linetype = "dashed") +
  facet_wrap(~ Time, ncol = 3) +
  labs(title = "Q-Q Plots: Total OxPLs by Time Point",
       x     = "Theoretical Quantiles",
       y     = "Sample Quantiles") +
  theme_bw() +
  theme(plot.title       = element_text(hjust = 0.5, face = "bold", size = 12),
        strip.background = element_rect(fill = "grey90"),
        strip.text       = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

print(qq_oxpl)

# 5.5 Q-Q plot: BODIPY C11
qq_bodipy <- ggplot(BODIPY_long_for_norm, aes(sample = Value)) +
  stat_qq(color = "#228B22", size = 2) +
  stat_qq_line(color = "black", linetype = "dashed") +
  facet_wrap(~ Time, ncol = 3) +
  labs(title = "Q-Q Plots: BODIPY C11 by Time Point",
       x     = "Theoretical Quantiles",
       y     = "Sample Quantiles") +
  theme_bw() +
  theme(plot.title       = element_text(hjust = 0.5, face = "bold", size = 12),
        strip.background = element_rect(fill = "grey90"),
        strip.text       = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

print(qq_bodipy)

# ============================================================
# 6. Statistical testing: each time point vs 0 min # -----
# ============================================================

baseline_values <- data_total %>%
  filter(Time == "0 min") %>%
  pull(Total_Intensity)

t_test_results <- data.frame()

for (tp in time_points) {
  current_values <- data_total %>%
    filter(Time == tp) %>%
    pull(Total_Intensity)
  
  test_result <- t.test(current_values, baseline_values,
                        paired    = FALSE,
                        var.equal = FALSE)
  
  t_test_results <- rbind(t_test_results, data.frame(
    Time           = tp,
    Mean_Baseline  = mean(baseline_values),
    Mean_TimePoint = mean(current_values),
    t_statistic    = test_result$statistic,
    df             = test_result$parameter,
    p_value        = test_result$p.value,
    CI_lower       = test_result$conf.int[1],
    CI_upper       = test_result$conf.int[2]
  ))
}

t_test_results$p_adj_BH    <- p.adjust(t_test_results$p_value, method = "BH")
t_test_results$Significance <- cut(t_test_results$p_adj_BH,
                                   breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                   labels = c("***", "**", "*", "ns"))

print("Statistical Comparison vs 0 min Baseline:")
print(t_test_results)

# ============================================================
# 7. BODIPY C11: Z-score scaling and summary # -----
# ============================================================

BODIPY_long <- BODIPY %>%
  pivot_longer(cols     = -Replicate,
               names_to  = "Time_raw",
               values_to = "Value") %>%
  filter(!is.na(Value)) %>%
  mutate(Time = gsub("X",    "", Time_raw),
         Time = gsub("\\.", " ", Time),
         Time = factor(Time, levels = time_order))

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

# ============================================================
# 8. Pearson correlation: Total OxPLs vs BODIPY # -----
# ============================================================

correlation_data <- data_total_summary %>%
  select(Time, OxPL_Mean = Mean_Total_Scaled) %>%
  left_join(BODIPY_stats %>% select(Time, BODIPY_Mean = Mean_Scaled), by = "Time")

pearson_test <- cor.test(correlation_data$OxPL_Mean,
                         correlation_data$BODIPY_Mean,
                         method = "pearson")

correlation_result <- data.frame(
  Comparison = "Total OxPLs vs BODIPY",
  Pearson_r  = pearson_test$estimate,
  Pearson_p  = pearson_test$p.value
) %>%
  mutate(
    Pearson_sig = case_when(
      Pearson_p < 0.001 ~ "***",
      Pearson_p < 0.01  ~ "**",
      Pearson_p < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

print("Pearson Correlation Analysis:")
print(correlation_result)

legend_labels <- c(
  sprintf("Overall total OxPLs (*R* = %.3f, *P* = %.3f)",
          correlation_result$Pearson_r,
          correlation_result$Pearson_p),
  "Peroxidative-positive cells by BODIPY C11"
)

# ============================================================
# 9. Build bracket and significance annotation data # -----
# ============================================================

last_time <- "270 min"

endpoint_values <- bind_rows(
  data_total_summary %>%
    filter(Time == last_time) %>%
    select(y_end = Mean_Total_Scaled) %>%
    mutate(Curve = "OxPLs"),
  BODIPY_stats %>%
    filter(Time == last_time) %>%
    select(y_end = Mean_Scaled) %>%
    mutate(Curve = "BODIPY")
)

bracket_data <- data.frame(
  y_oxpl   = endpoint_values %>% filter(Curve == "OxPLs")  %>% pull(y_end),
  y_bodipy = endpoint_values %>% filter(Curve == "BODIPY") %>% pull(y_end)
) %>%
  mutate(
    y_min        = pmin(y_oxpl, y_bodipy),
    y_max        = pmax(y_oxpl, y_bodipy),
    y_mid        = (y_min + y_max) / 2,
    Significance = correlation_result$Pearson_sig,
    text_size    = case_when(Significance == "ns" ~ 3.5, TRUE ~ 6)
  )

x_bracket_base    <- 6.5
bracket_offset    <- 0
bracket_linewidth <- 0.3

bracket_segments <- data.frame(
  x    = c(x_bracket_base + bracket_offset,
           x_bracket_base + bracket_offset - 0.2,
           x_bracket_base + bracket_offset - 0.2),
  xend = c(x_bracket_base + bracket_offset,
           x_bracket_base + bracket_offset,
           x_bracket_base + bracket_offset),
  y    = c(bracket_data$y_min, bracket_data$y_max, bracket_data$y_min),
  yend = c(bracket_data$y_max, bracket_data$y_max, bracket_data$y_min),
  segment_type = c("vertical", "top", "bottom")
)

bracket_data <- bracket_data %>%
  mutate(x_label = x_bracket_base + bracket_offset + 0.1)

annotation_data <- data_total_summary %>%
  left_join(t_test_results %>% select(Time, p_adj_BH, Significance), by = "Time") %>%
  filter(Time != "0 min") %>%
  mutate(y_position  = Mean_Total_Scaled + SEM_Total_Scaled + 0.15,
         label_color = "#7E1FF2")

# ============================================================
# 10. Overlay plot: Z-scaled Total OxPLs + BODIPY C11 # -----
# ============================================================

p_total_annotated <- ggplot(data_total_summary,
                            aes(x = Time, y = Mean_Total_Scaled, group = 1)) +
  geom_ribbon(aes(ymin = Mean_Total_Scaled - SEM_Total_Scaled,
                  ymax = Mean_Total_Scaled + SEM_Total_Scaled,
                  fill = "Total OxPLs"),
              alpha = 0.2) +
  geom_line(aes(color = "Total OxPLs"), linewidth = 1.2) +
  geom_point(aes(color = "Total OxPLs", shape = "Total OxPLs"), size = 4) +
  geom_ribbon(data = BODIPY_stats,
              aes(x = Time, y = Mean_Scaled, group = 1,
                  ymin = Mean_Scaled - SEM,
                  ymax = Mean_Scaled + SEM,
                  fill = "BODIPY"),
              inherit.aes = FALSE, alpha = 0.2) +
  geom_line(data = BODIPY_stats,
            aes(x = Time, y = Mean_Scaled, group = 1, color = "BODIPY"),
            inherit.aes = FALSE, linewidth = 1.2) +
  geom_point(data = BODIPY_stats,
             aes(x = Time, y = Mean_Scaled, color = "BODIPY", shape = "BODIPY"),
             inherit.aes = FALSE, size = 4) +
  geom_text(data = annotation_data,
            aes(x = Time, y = y_position, label = Significance),
            size = 4, vjust = 0, color = "#7E1FF2") +
  geom_segment(data = bracket_segments,
               aes(x = x, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE, color = "black", linewidth = bracket_linewidth) +
  geom_text(data = bracket_data,
            aes(x = x_label, y = y_mid, label = Significance, size = text_size),
            inherit.aes = FALSE, color = "black",
            hjust = 0, vjust = 0.5, fontface = "plain") +
  scale_size_identity() +
  scale_x_discrete(labels = setNames(time_labels_clean, time_order)) +
  scale_y_continuous(limits = c(-1.7, 2.7), expand = expansion(mult = c(0, 0.05))) +
  scale_color_manual(values = c("BODIPY" = "#228B22", "Total OxPLs" = "#7E1FF2"),
                     labels = c("BODIPY"      = legend_labels[2],
                                "Total OxPLs" = legend_labels[1]),
                     breaks = c("Total OxPLs", "BODIPY"),
                     name = "") +
  scale_fill_manual(values = c("BODIPY" = "#228B22", "Total OxPLs" = "#7E1FF2"),
                    labels = c("BODIPY"      = legend_labels[2],
                               "Total OxPLs" = legend_labels[1]),
                    breaks = c("Total OxPLs", "BODIPY"),
                    name = "") +
  scale_shape_manual(values = c("BODIPY" = 16, "Total OxPLs" = 17),
                     labels = c("BODIPY"      = legend_labels[2],
                                "Total OxPLs" = legend_labels[1]),
                     breaks = c("Total OxPLs", "BODIPY"),
                     name = "") +
  guides(fill  = guide_legend(override.aes = list(size = c(3.5, 3.5)),
                              keyheight = unit(1.3, "lines")),
         shape = guide_legend(override.aes = list(size = c(3.5, 3.5)),
                              keyheight = unit(1.3, "lines")),
         color = guide_legend(override.aes = list(size = c(3.5, 3.5)),
                              keyheight = unit(1.3, "lines"))) +
  coord_cartesian(clip = "off", xlim = c(0.8, 6.5)) +
  theme_bw() +
  theme(
    axis.text.x      = element_text(angle = 0, hjust = 0.5, size = 13),
    axis.text.y      = element_text(size = 13),
    axis.title       = element_text(size = 14, face = "plain"),
    axis.title.x     = element_text(size = 14, face = "plain", vjust = -1.5),
    plot.title       = element_markdown(size = 16, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", linewidth = 0.3),
    legend.position  = "right",
    legend.text      = element_markdown(size = 11),
    legend.title     = element_text(size = 14, face = "bold"),
    legend.key.size  = unit(1.3, "lines"),
    plot.margin      = unit(c(5, 5, 5, 5), "mm")
  ) +
  labs(x = "ML210 treatment time (min)",
       y = "Z score")

print(p_total_annotated)

# ============================================================
# 11. Export figures # -----
# ============================================================

pdf("./Figure/overlay_of_total_OxPLs_and_Bodipy.pdf", width = 750 / 96, height = 320 / 96)
print(p_total_annotated)
dev.off()

pdf("./Figure/Normality_QQ_TotalOxPLs.pdf", width = 600 / 96, height = 400 / 96)
print(qq_oxpl)
dev.off()

pdf("./Figure/Normality_QQ_BODIPY.pdf", width = 600 / 96, height = 400 / 96)
print(qq_bodipy)
dev.off()