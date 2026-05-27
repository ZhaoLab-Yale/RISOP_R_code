rm(list = ls())
library(dplyr)
library(ggplot2)
library(signal)

# =====================================================
# 1. Fenton ----
# =====================================================
PS_Fenton <- read.csv('./Input/20260221_Neg_FentonQC_PRM_1_PS_EIC.csv', stringsAsFactors = F, skip = 3, header = T)

start_time <- 6.4
end_time <- 8.8

df_plot <- data.frame(
  Time = as.numeric(PS_Fenton$Time.min.),
  Intensity = as.numeric(PS_Fenton$Intensity)
)

df_plot <- df_plot[!is.na(df_plot$Time) & !is.na(df_plot$Intensity), ]
df_plot <- df_plot[order(df_plot$Time), ]
df_plot <- df_plot[df_plot$Time >= start_time & df_plot$Time <= end_time, ]

# Subtract mean intensity of a flat baseline region to correct for background signal
baseline_region <- df_plot[df_plot$Time >= 6.5 & df_plot$Time <= 7, ]
baseline <- mean(baseline_region$Intensity, na.rm = TRUE)
df_plot$Intensity <- df_plot$Intensity - baseline

df_plot$Intensity[df_plot$Intensity < 0] <- 0

# LOESS smoothing — handles asymmetric peaks and tails well
loess_model <- loess(Intensity ~ Time, data = df_plot, span = 0.18, degree = 0)
df_plot$SmoothIntensity <- predict(loess_model, df_plot$Time)

df_plot$SmoothIntensity[df_plot$SmoothIntensity < 0] <- 0

a <- ggplot(df_plot, aes(x = Time)) +
  #geom_line(aes(y = Intensity), color = "grey70", linewidth = 0.5, alpha = 0.5) +
  geom_line(aes(y = SmoothIntensity), color = "#E05A51", linewidth = 1) +
  scale_y_continuous(limits = c(0, 50000), expand = c(0, 0), 
                     breaks = seq(0, 50000, by = 10000)) +
  scale_x_continuous(limits = c(start_time, end_time), expand = c(0.01, 0)) +
  labs(
    x = "Retention time (min)",
    y = "Intensity"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 0, unit = "pt")
  )

ggsave("./Figure/PS_Fenton_EIC.pdf", plot = a, width = 300/96, height = 250/96, device = "pdf")

# =====================================================
# 2. H2O2 ----
# =====================================================
PS_H2O2 <- read.csv('./Input/20260221_Neg_H2O2QC_PRM_1_PS_EIC.csv', stringsAsFactors = F, skip = 3, header = T)

start_time <- 6.4
end_time <- 8.8

df_plot <- data.frame(
  Time = as.numeric(PS_H2O2$Time.min.),
  Intensity = as.numeric(PS_H2O2$Intensity)
)

df_plot <- df_plot[!is.na(df_plot$Time) & !is.na(df_plot$Intensity), ]
df_plot <- df_plot[order(df_plot$Time), ]
df_plot <- df_plot[df_plot$Time >= start_time & df_plot$Time <= end_time, ]

# Subtract mean intensity of a flat baseline region to correct for background signal
baseline_region <- df_plot[df_plot$Time >= 6.5 & df_plot$Time <= 7, ]
baseline <- mean(baseline_region$Intensity, na.rm = TRUE)
df_plot$Intensity <- df_plot$Intensity - baseline

df_plot$Intensity[df_plot$Intensity < 0] <- 0

# LOESS smoothing — handles asymmetric peaks and tails well
loess_model <- loess(Intensity ~ Time, data = df_plot, span = 0.18, degree = 0)
df_plot$SmoothIntensity <- predict(loess_model, df_plot$Time)

df_plot$SmoothIntensity[df_plot$SmoothIntensity < 0] <- 0

b <- ggplot(df_plot, aes(x = Time)) +
  #geom_line(aes(y = Intensity), color = "grey70", linewidth = 0.5, alpha = 0.5) +
  geom_line(aes(y = SmoothIntensity), color = "#377EB8", linewidth = 1) +
  scale_y_continuous(limits = c(0, 50000), expand = c(0, 0), 
                     breaks = seq(0, 50000, by = 10000)) +
  scale_x_continuous(limits = c(start_time, end_time), expand = c(0.01, 0)) +
  labs(
    x = "Retention time (min)",
    y = "Intensity"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 0, unit = "pt")
  )

ggsave("./Figure/PS_H2O2_EIC.pdf", plot = b, width = 300/96, height = 250/96, device = "pdf")