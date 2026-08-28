rm(list = ls())
library(dplyr)
library(ggplot2)

# ============================================================
# 1. Load data ----
# ============================================================
PC     <- read.csv('./Input/20260221_Neg_FentonQC_PRM_1_PC_MS2.csv', stringsAsFactors = F, skip = 7, header = T)
PS     <- read.csv('./Input/20260221_Neg_H2O2QC_PRM_1_PS_MS2.csv', stringsAsFactors = F, skip = 7, header = T)
PC_ref <- read.csv('./Input/PC_ref.csv', stringsAsFactors = F, skip = 7, header = T)
PS_ref <- read.csv('./Input/PS_ref.csv', stringsAsFactors = F, skip = 7, header = T)

# ============================================================
# 2. Process MS2 data ----
# ============================================================
# ??? Change to PC/PS pair:
top_raw    <- PS
bottom_raw <- PS_ref

process_ms2 <- function(raw, sign = 1) {
  raw %>%
    mutate(
      Mass      = as.numeric(Mass),
      Intensity = as.numeric(Intensity)
    ) %>%
    filter(!is.na(Mass), !is.na(Intensity)) %>%
    arrange(Mass) %>%
    mutate(RelIntensity = sign * 100 * Intensity / max(Intensity, na.rm = TRUE))
}

df_top    <- process_ms2(top_raw,    sign = 1)  %>% mutate(Source = "Representative")
df_bottom <- process_ms2(bottom_raw, sign = -1) %>% mutate(Source = "Reference")

df <- bind_rows(df_top, df_bottom)

# ============================================================
# 3. Plot mirror MS2 spectrum ----
# ============================================================
x_range <- range(df$Mass)

p <- ggplot(df, aes(x = Mass, y = RelIntensity)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  geom_segment(aes(xend = Mass, y = 0, yend = RelIntensity),
               linewidth = 0.6, color = "black") +
  labs(
    x = "m/z",
    y = "Relative abundance"
  ) +
  scale_x_continuous(
    breaks = seq(floor(x_range[1] / 50) * 50, ceiling(x_range[2] / 50) * 50, by = 100),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(-100, 100),
    breaks = seq(-100, 100, by = 50),
    labels = abs(seq(-100, 100, by = 50)),
    expand = c(0, 0)
  ) +
  theme_classic() +
  theme(
    panel.grid  = element_blank(),
    panel.border = element_blank(),
    axis.line   = element_line(color = "#787878", linewidth = 0.4),
    axis.ticks  = element_line(color = "#787878", linewidth = 0.4),
    text        = element_text(family = "Arial"),
    axis.title  = element_text(size = 14, face = "bold"),
    axis.text   = element_text(size = 10, color = "black"),
    plot.title  = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

p

ggsave("./Output/PS_spectra.png", plot = p, width = 6, height = 3.75, dpi = 300, bg = "white")