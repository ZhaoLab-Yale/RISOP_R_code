rm(list = ls())
library(tidyverse)
library(scales)
library(openxlsx)
library(writexl)
library(readxl)
setwd(rstudioapi::getActiveProject())

# ============================================================
# 1. Load data # -----
# ============================================================

Ox_ML210 <- read.xlsx("./Output/Table S5. Normalized intensities of OxPLs and molar concentrations of non-oxidized lipids identified in ML210-treated cells.xlsx", sheet = "OxPLs (filtered)") %>%
  select(Lipid, Identification.source)

# ============================================================
# 2. Compute counts # -----
# ============================================================

# 2.1 Count lipids per identification source combination
counts <- Ox_ML210 %>%
  summarise(
    n_samples      = sum(str_detect(Identification.source, "ML210")),
    n_fenton_bar2  = sum(str_detect(Identification.source, "Fenton") & !str_detect(Identification.source, "ML210")),
    n_h2o2_bar3    = sum(str_detect(Identification.source, "H2O2")   & !str_detect(Identification.source, "ML210")),
    n_fenton_bar4  = sum(str_detect(Identification.source, "Fenton") & !str_detect(Identification.source, "ML210") & !str_detect(Identification.source, "H2O2")),
    n_h2o2_bar4    = sum(str_detect(Identification.source, "H2O2")   & !str_detect(Identification.source, "ML210") & !str_detect(Identification.source, "Fenton")),
    n_overlap_bar4 = sum(str_detect(Identification.source, "Fenton") &  str_detect(Identification.source, "H2O2")  & !str_detect(Identification.source, "ML210"))
  )

# 2.2 Assemble stacked bar data frame
df <- data.frame(
  Group = c(
    rep("ML210",             4),
    rep("ML210+Fenton",      4),
    rep("ML210+H2O2",        4),
    rep("ML210+Fenton+H2O2", 4)
  ),
  FillCat = rep(c("Samples", "Overlap", "Fenton", "H2O2"), 4),
  n = c(
    counts$n_samples, 0, 0, 0,
    counts$n_samples, 0, counts$n_fenton_bar2, 0,
    counts$n_samples, 0, 0, counts$n_h2o2_bar3,
    counts$n_samples, counts$n_overlap_bar4, counts$n_fenton_bar4, counts$n_h2o2_bar4
  )
) %>%
  mutate(
    Group   = factor(Group,   levels = c("ML210+Fenton+H2O2", "ML210+H2O2", "ML210+Fenton", "ML210")),
    FillCat = factor(FillCat, levels = c("Samples", "Fenton", "H2O2", "Overlap"))
  )

# ============================================================
# 3. Define colors # -----
# ============================================================

fill_cols <- c(
  "Samples" = alpha("#7E1FF2", 0.8),
  "Fenton"  = alpha("#E05A51", 1),
  "H2O2"    = alpha("#377EB8", 0.9),
  "Overlap" = alpha("#6B1C4D", 0.8)
)

# ============================================================
# 4. Plot and export # -----
# ============================================================

# 4.1 Build stacked bar plot
p_final <- ggplot(df, aes(x = Group, y = n, fill = FillCat)) +
  geom_col(
    width     = 0.75,
    color     = NA,
    linewidth = 0,
    position  = position_stack(reverse = TRUE)
  ) +
  scale_fill_manual(
    values = fill_cols,
    labels = c(
      "Samples" = expression("ML210-treated samples"),
      "Fenton"  = expression("Additional OxPLs identified from Fenton references"),
      "H2O2"    = expression("Additional OxPLs identified from" ~ H[2]*O[2] * "-only references"),
      "Overlap" = expression("Additional OxPLs identified from both Fenton and" ~ H[2]*O[2] * "-only references")
    )
  ) +
  scale_x_discrete(
    labels = c(
      "ML210"             = "Samples alone",
      "ML210+Fenton"      = "Samples + Fenton references",
      "ML210+H2O2"        = expression(Samples ~ "+" ~ H[2]*O[2]*"-only" ~ references),
      "ML210+Fenton+H2O2" = "Samples + both references"
    )
  ) +
  scale_y_continuous(
    breaks = seq(0, 25, by = 5),
    expand = expansion(mult = c(0, 0.05))
  ) +
  coord_flip() +
  labs(x = NULL, y = "Count", fill = NULL) +
  guides(
    fill = guide_legend(
      keywidth  = unit(0.6, "cm"),
      keyheight = unit(0.55, "cm"),
      byrow     = TRUE
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.line    = element_line(linewidth = 0.5),
    axis.ticks   = element_line(linewidth = 0.5),
    axis.text.y  = element_text(size = 21, hjust = 1),
    axis.text.x  = element_text(size = 20),
    axis.title.x = element_text(face = "bold", size = 20),
    
    plot.title.position = "plot",
    plot.title          = element_text(size = 24, face = "bold", hjust = 0.3),
    
    legend.position      = "right",
    legend.title         = element_blank(),
    legend.text          = element_text(size = 19, lineheight = 1.1),
    legend.key.spacing.y = unit(0.4, "cm"),
    plot.margin          = margin(10, 5, 10, 5)
  )

p_final

# 4.2 Export as SVG
ggsave(
  "./Figure/ML210_column.svg",
  plot   = p_final,
  width  = 1700/96, height = 415/96, units = "in",
  device = svglite::svglite,
  bg     = "transparent"
)