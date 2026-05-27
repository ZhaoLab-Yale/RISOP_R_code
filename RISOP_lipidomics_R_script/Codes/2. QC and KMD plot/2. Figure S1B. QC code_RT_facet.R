rm(list = ls())
library(tidyverse)
library(ggplot2)
library(readxl)
setwd(rstudioapi::getActiveProject())

# ==============================================================================
# 1. Load input data -----
# ==============================================================================
rt_data <- read_xlsx('./Input/KP4 retention time.xlsx', sheet = 1)
colnames(rt_data)[1] <- "Compound"

# ==============================================================================
# 2. Data preparation -----
# ==============================================================================

# 2.1 Reshape to long format and set QC injection order as ordered factor
rt_long <- rt_data %>%
  pivot_longer(cols = -Compound, names_to = "QC", values_to = "RetentionTime") %>%
  mutate(QC = factor(QC, levels = c("#1", "#4", "#14", "#26", "#38")))

# 2.2 Calculate mean RT, CV, and ±0.08 / ±0.05 min reference band per compound
rt_stats <- rt_long %>%
  group_by(Compound) %>%
  summarise(
    mean_RT = mean(RetentionTime, na.rm = TRUE),
    cv_pct  = (sd(RetentionTime, na.rm = TRUE) / mean(RetentionTime, na.rm = TRUE)) * 100,
    y_min   = mean_RT - 0.08,
    y_max   = mean_RT + 0.08,
    y_min10 = mean_RT - 0.05,  # Inner ±0.05 min tolerance band
    y_max10 = mean_RT + 0.05,
    .groups = "drop"
  )

# 2.3 Dummy data to force free_y axis range to ±0.08 min around mean
dummy_data <- rt_stats %>%
  select(Compound, y_min, y_max) %>%
  pivot_longer(cols = c(y_min, y_max), values_to = "RetentionTime") %>%
  mutate(QC = factor("#1", levels = levels(rt_long$QC)))

# 2.4 CV annotation labels per compound
cv_labels <- rt_stats %>%
  mutate(label = paste0("CV: ", round(cv_pct, 2), "%"))

# ==============================================================================
# 3. Retention time QC plot -----
# ==============================================================================
p <- ggplot(rt_long, aes(x = QC, y = RetentionTime, group = Compound)) +
  # Force y-axis range via invisible dummy points
  geom_blank(data = dummy_data) +
  # ±0.05 min tolerance band reference lines
  geom_hline(data = rt_stats, aes(yintercept = y_max10),
             linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_hline(data = rt_stats, aes(yintercept = y_min10),
             linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  # Retention time trajectory per compound
  geom_line(linewidth = 0.7, colour = "#1F4E79") +
  # Individual QC injection points
  geom_point(size = 2.5,
             shape  = 21,
             colour = "#4C78A8",
             fill   = "#4C78A8",
             alpha  = 0.9,
             stroke = 0.6) +
  # CV annotation in top-right corner of each facet
  geom_text(data = cv_labels,
            aes(x = Inf, y = Inf, label = label),
            hjust = 1.25, vjust = 10,
            size = 3.2, colour = "black", inherit.aes = FALSE) +
  facet_wrap(~ Compound, scales = "free_y", ncol = 4) +
  labs(
    x = "Injection order of pooled sample",
    y = "Retention time (min)"
  ) +
  theme_bw() +
  theme(
    axis.text.x      = element_text(size = 11, angle = 0, hjust = 0.5),
    axis.text.y      = element_text(size = 11),
    axis.title.x     = element_text(size = 15.23, face = "bold", margin = margin(t = 12, b = 0)),
    axis.title.y     = element_text(size = 15.23, face = "bold", margin = margin(r = 16)),
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text       = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.y  = unit(0.8, "cm"),
    plot.margin      = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")
  )

print(p)

# ==============================================================================
# 4. Export -----
# ==============================================================================
pdf("./Figure/QC_retention_time.pdf", width = 1000/96, height = 575/96)
print(p)
dev.off()