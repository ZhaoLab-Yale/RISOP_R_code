rm(list = ls())
library(dplyr)
library(ggplot2)

# ============================================================
# 1. Load data # -----
# ============================================================

PC <- read.csv('./Input/20260221_Neg_FentonQC_PRM_1_PC_MS2.csv', stringsAsFactors = F, skip = 7, header = T)
PS <- read.csv('./Input/20260221_Neg_H2O2QC_PRM_1_PS_MS2.csv', stringsAsFactors = F, skip = 7, header = T)

# ============================================================
# 2. Process MS2 data # -----
# ============================================================

# ??? Change to PC or PS:
df <- PS %>%
  mutate(
    Mass      = as.numeric(Mass),
    Intensity = as.numeric(Intensity)
  ) %>%
  filter(!is.na(Mass), !is.na(Intensity)) %>%
  arrange(Mass) %>%
  mutate(RelIntensity = 100 * Intensity / max(Intensity, na.rm = TRUE))

# ============================================================
# 3. Plot MS2 spectrum # -----
# ============================================================

ggplot(df, aes(x = Mass, y = RelIntensity)) +
  geom_segment(aes(xend = Mass, y = 0, yend = RelIntensity),
               linewidth = 0.6,
               color = "black") +
  labs(
    x = "m/z",
    y = "Relative Intensity"
  ) +
  scale_x_continuous(
    breaks = seq(
      floor(min(df$Mass) / 50) * 50,
      ceiling(max(df$Mass) / 50) * 50,
      by = 100
    ),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  theme_classic() +
  theme(
    panel.grid  = element_blank(),
    panel.border = element_blank(),
    # Axes
    axis.line   = element_line(color = "#787878", linewidth = 0.4),
    axis.ticks  = element_line(color = "#787878", linewidth = 0.4),
    # Fonts
    text        = element_text(family = "Arial"),
    axis.title  = element_text(size = 14, face = "bold"),
    axis.text   = element_text(size = 10, color = "black"),
    plot.title  = element_text(size = 16, face = "bold", hjust = 0.5)
  )