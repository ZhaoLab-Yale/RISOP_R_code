rm(list = ls())
library(tidyverse)
library(readxl)
library(ggplot2)
library(drc)
library(patchwork)
library(cowplot)
library(grid)
setwd(rstudioapi::getActiveProject())

# ============================================================
# 1. Shared plot settings # -----
# ============================================================

x_limits <- c(-20, 300)

common_x <- scale_x_continuous(
  breaks = seq(0, 300, by = 60),
  limits = x_limits,
  expand = c(0, 0)
)

common_theme <- theme_classic() %+replace%
  theme(
    axis.text.x       = element_text(size = 13, face = "bold", color = "black"),
    axis.text.y       = element_text(size = 13, face = "bold", color = "black"),
    axis.line         = element_line(color = "black", linewidth = 0.8),
    axis.ticks        = element_line(color = "black", linewidth = 0.8),
    axis.ticks.length = unit(2.5, "mm"),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank()
  )

# ============================================================
# 2. Broken-axis helper functions # -----
# ============================================================

# 2.1 Build two parallel slash grobs to mark axis break
make_parallel_slashes_fixed <- function(angle_deg, length_mm,
                                        gap_mm, lwd_pt,
                                        x_offset_mm, y_npc) {
  angle_rad <- angle_deg * pi / 180
  dx <- unit(length_mm * cos(angle_rad), "mm")
  dy <- unit(length_mm * sin(angle_rad), "mm")
  
  x0   <- unit(0, "npc") - unit(x_offset_mm, "mm")
  y0_1 <- unit(y_npc, "npc")
  y0_2 <- unit(y_npc, "npc") + unit(gap_mm, "mm")
  
  make_one <- function(y_start, col, lwd) {
    segmentsGrob(
      x0 = x0,
      y0 = y_start,
      x1 = x0 + dx,
      y1 = y_start + dy,
      gp = gpar(col = col, lwd = lwd, lineend = "butt")
    )
  }
  
  bg1 <- make_one(y0_1, "white", lwd_pt + 2)
  fg1 <- make_one(y0_1, "black", lwd_pt)
  bg2 <- make_one(y0_2, "white", lwd_pt + 2)
  fg2 <- make_one(y0_2, "black", lwd_pt)
  
  grobTree(bg1, fg1, bg2, fg2)
}

# 2.2 Find lower panel viewport name dynamically
find_lower_panel_viewport <- function() {
  vp_str    <- capture.output(print(current.vpTree(all = TRUE)))
  all_names <- regmatches(vp_str,
                          gregexpr("(?<=viewport\\[)[^\\]]+", vp_str, perl = TRUE))
  all_names <- unique(unlist(all_names))
  
  patchwork_panels <- all_names[grepl("^panel-[0-9]+\\.", all_names)]
  
  if (length(patchwork_panels) > 0) {
    extract_idx <- function(nm) {
      m <- regmatches(nm, regexpr("(?<=panel-)[0-9]+", nm, perl = TRUE))
      if (length(m) == 0) return(0L) else return(as.integer(m))
    }
    idxs     <- sapply(patchwork_panels, extract_idx)
    lower_vp <- patchwork_panels[which.max(idxs)]
    return(lower_vp)
  }
  
  panel_names <- all_names[grepl("^panel", all_names)]
  if (length(panel_names) == 0)
    stop("No panel viewports found. All names:\n", paste(all_names, collapse = "\n"))
  
  extract_first_int <- function(nm) {
    m <- regmatches(nm, regexpr("[0-9]+", nm))
    if (length(m) == 0) return(0L) else return(as.integer(m))
  }
  rows     <- sapply(panel_names, extract_first_int)
  lower_vp <- panel_names[which.max(rows)]
  return(lower_vp)
}

# 2.3 Draw final figure with slashes anchored to the lower panel viewport
draw_final_with_slashes <- function(plot_obj,
                                    angle_deg, length_mm,
                                    gap_mm, lwd_pt,
                                    x_offset_mm = 3.0,
                                    slash_y_npc = 1.0) {
  gt <- as_gtable(plot_obj)
  
  grid.newpage()
  pushViewport(viewport(width = 1, height = 1))
  grid.draw(gt)
  
  lower_vp <- find_lower_panel_viewport()
  message("Using viewport: ", lower_vp)
  
  downViewport(lower_vp)
  
  slash_grob <- make_parallel_slashes_fixed(
    angle_deg   = angle_deg,
    length_mm   = length_mm,
    gap_mm      = gap_mm,
    lwd_pt      = lwd_pt,
    x_offset_mm = x_offset_mm,
    y_npc       = slash_y_npc
  )
  
  grid.draw(slash_grob)
  upViewport(0)
}

# ============================================================
# 3. Cell viability over time (Sheet 1) # -----
# ============================================================

# 3.1 Adjust parameters
ylab_x <- 0.05
ylab_y <- 0.6

slash_angle_deg   <- 30
slash_length      <- 5
slash_gap_mm      <- 2.0
slash_lwd_pt      <- 2
slash_x_offset_mm <- 2.3
slash_y_npc       <- 0.96

xlab_margin_t  <- 10
xlab_margin_b  <- 5
plot_bottom_mm <- 0

# 3.2 Load and reshape data
df_raw <- read_excel("./Input/Cell viability.xlsx", sheet = 1)
colnames(df_raw) <- c("Time", "Rep1", "Rep2", "Rep3")

df_long <- df_raw %>%
  pivot_longer(
    cols      = c(Rep1, Rep2, Rep3),
    names_to  = "Replicate",
    values_to = "Viability"
  )

df_summary <- df_long %>%
  group_by(Time) %>%
  summarise(
    Mean_Viability = mean(Viability, na.rm = TRUE),
    SEM            = sd(Viability,   na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# 3.3 Upper panel (y: 97–100.3)
p_upper <- ggplot(df_summary, aes(x = Time, y = Mean_Viability)) +
  geom_line(color = "#006400", linewidth = 1.2) +
  geom_errorbar(
    aes(ymin = Mean_Viability - SEM, ymax = Mean_Viability + SEM),
    width = 11, color = "black", linewidth = 0.5
  ) +
  geom_point(color = "black", shape = 16, size = 1.5) +
  common_x +
  scale_y_continuous(
    limits = c(97, 100.3),
    breaks = c(98, 100),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  common_theme +
  theme(
    axis.line.x.bottom  = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.text.x.bottom  = element_blank(),
    axis.title.x        = element_blank(),
    axis.title.y        = element_blank(),
    plot.margin         = unit(c(2, 10, 1, 10), "mm")
  )

# 3.4 Lower panel (y: 0–88)
p_lower <- ggplot(df_summary, aes(x = Time, y = Mean_Viability)) +
  geom_line(color = "#006400", linewidth = 1.2) +
  geom_errorbar(
    aes(ymin = Mean_Viability - SEM, ymax = Mean_Viability + SEM),
    width = 8, color = "black", linewidth = 0.8
  ) +
  geom_point(color = "black", shape = 16, size = 2.2) +
  common_x +
  scale_y_continuous(
    limits = c(0, 88),
    breaks = c(0, 25, 50, 75),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  common_theme +
  theme(
    axis.title.x     = element_text(size = 14.5, face = "bold",
                                    margin = margin(t = xlab_margin_t,
                                                    b = xlab_margin_b)),
    axis.title.y     = element_blank(),
    axis.line.x.top  = element_blank(),
    axis.ticks.x.top = element_blank(),
    plot.margin      = unit(c(1, 10, plot_bottom_mm, 10), "mm")
  ) +
  labs(x = "ML210 treatment time (min)")

# 3.5 Stack panels and add shared y-axis label
p_stacked <- (p_upper / p_lower) +
  plot_layout(heights = c(1.3, 2.7))

p_with_ylab <- ggdraw(p_stacked) +
  draw_label(
    "Cell viability (%)",
    x        = ylab_x,
    y        = ylab_y,
    angle    = 90,
    fontface = "bold",
    size     = 14.5,
    hjust    = 0.5
  )

# 3.6 Preview and export
draw_final_with_slashes(
  plot_obj    = p_with_ylab,
  angle_deg   = slash_angle_deg,
  length_mm   = slash_length,
  gap_mm      = slash_gap_mm,
  lwd_pt      = slash_lwd_pt,
  x_offset_mm = slash_x_offset_mm,
  slash_y_npc = slash_y_npc
)

quartz(type = "pdf", file = "./Figure/Cell_viability_time.pdf",
       width = 430 / 96, height = 280 / 96)
draw_final_with_slashes(
  plot_obj    = p_with_ylab,
  angle_deg   = slash_angle_deg,
  length_mm   = slash_length,
  gap_mm      = slash_gap_mm,
  lwd_pt      = slash_lwd_pt,
  x_offset_mm = slash_x_offset_mm,
  slash_y_npc = slash_y_npc
)
dev.off()

# ============================================================
# 4. ML210 dose-response (Sheet 2) # -----
# ============================================================

# 4.1 Load and reshape data
df_raw2 <- read_excel("./Input/Cell viability.xlsx", sheet = 2)

colnames(df_raw2) <- c(
  "Concentration",
  "W1", "W2", "W3", "W4", "W5", "W6",
  "F1", "F2", "F3", "F4", "F5", "F6"
)

df_raw2 <- df_raw2 %>%
  filter(Concentration > 0)

df_long2 <- df_raw2 %>%
  pivot_longer(
    cols      = -Concentration,
    names_to  = "Replicate",
    values_to = "Viability"
  ) %>%
  mutate(
    Group = ifelse(
      Replicate %in% c("W1", "W2", "W3", "W4", "W5", "W6"),
      "Without Fer-1",
      "With Fer-1"
    ),
    Group = factor(Group, levels = c("With Fer-1", "Without Fer-1"))
  )

df_summary2 <- df_long2 %>%
  group_by(Concentration, Group) %>%
  summarise(
    Mean_Viability = mean(Viability, na.rm = TRUE),
    SEM            = sd(Viability,   na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# 4.2 Fit dose-response curves
conc_seq <- exp(seq(
  log(min(df_raw2$Concentration)),
  log(max(df_raw2$Concentration)),
  length.out = 200
))

fit_without <- drm(
  Viability ~ Concentration,
  data = df_long2 %>% filter(Group == "Without Fer-1"),
  fct  = LL.4(names = c("Slope", "Lower", "Upper", "EC50"))
)

curve_without <- data.frame(
  Concentration = conc_seq,
  Viability     = predict(fit_without, newdata = data.frame(Concentration = conc_seq)),
  Group         = "Without Fer-1"
)

curve_with <- data.frame(
  Concentration = conc_seq,
  Viability     = 100,
  Group         = "With Fer-1"
)

curve_all <- bind_rows(curve_without, curve_with) %>%
  mutate(Group = factor(Group, levels = c("With Fer-1", "Without Fer-1")))

# 4.3 Build dose-response plot
p2 <- ggplot() +
  geom_line(
    data = curve_all,
    aes(x = Concentration, y = Viability, color = Group),
    linewidth = 1.2
  ) +
  geom_errorbar(
    data = df_summary2,
    aes(x = Concentration, ymin = Mean_Viability - SEM, ymax = Mean_Viability + SEM),
    width = 0.15, color = "black", linewidth = 0.6
  ) +
  geom_point(
    data = df_summary2,
    aes(x = Concentration, y = Mean_Viability),
    color = "black", shape = 16, size = 1.5
  ) +
  scale_x_log10(
    breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
    labels = c("0.001", "0.01", "0.1", "1", "10", "100"),
    limits = c(0.003, 30)
  ) +
  scale_y_continuous(
    limits = c(0, 130),
    breaks = seq(0, 125, by = 25),
    expand = c(0, 0)
  ) +
  scale_color_manual(
    values = c(
      "With Fer-1"    = "#7D5C47",
      "Without Fer-1" = "#C71585"
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(
    axis.text.x        = element_text(size = 13, face = "bold", color = "black"),
    axis.text.y        = element_text(size = 13, face = "bold", color = "black"),
    axis.title.x       = element_text(size = 14, face = "bold", margin = margin(t = 8)),
    axis.title.y       = element_text(size = 14, face = "bold", margin = margin(r = 8)),
    axis.line          = element_line(color = "black", linewidth = 0.8),
    axis.ticks         = element_line(color = "black", linewidth = 0.8),
    axis.ticks.length  = unit(2.5, "mm"),
    legend.title       = element_blank(),
    legend.text        = element_text(size = 11, face = "bold"),
    legend.position    = "right",
    legend.box.spacing = unit(2, "mm"),
    legend.margin      = margin(0, 0, 0, -8),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.margin        = unit(c(5, 10, 5, 5), "mm")
  ) +
  labs(
    x = "ML210 (\u03bcM)",
    y = "Cell viability (%)"
  )

# 4.4 Preview and export
print(p2)

quartz(type = "pdf", file = "./Figure/Cell_viability_dose.pdf",
       width = 550 / 96, height = 287 / 96)
print(p2)
dev.off()