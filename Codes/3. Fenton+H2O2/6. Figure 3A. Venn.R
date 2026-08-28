rm(list = ls())
library(readxl)
library(VennDiagram)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
setwd(rstudioapi::getActiveProject())

# ============================================================
# 1. Load data # -----
# ============================================================

load(file = paste0("./Output/cmbd_lipids_Fenton.Rdata"))
load(file = paste0("./Output/cmbd_lipids_H2O2.Rdata"))

# ============================================================
# 2. Process lipid lists # -----
# ============================================================

cmbd_lipids_Fenton <- cmbd_lipids_Fenton %>%
  filter(grepl("<", Metabolite.curated, fixed = TRUE)) %>%
  select(Metabolite.curated)

cmbd_lipids_H2O2 <- cmbd_lipids_H2O2 %>%
  filter(grepl("<", Metabolite.curated, fixed = TRUE)) %>%
  select(Metabolite.curated)

# ============================================================
# 3. Plot and export Venn diagram # -----
# ============================================================

# 3.1 Build sets
fenton_set <- cmbd_lipids_Fenton$Metabolite.curated
h2o2_set   <- cmbd_lipids_H2O2$Metabolite.curated

# 3.2 Draw Venn diagram
venn.plot <- venn.diagram(
  x = list(
    Fenton = fenton_set,
    H2O2   = h2o2_set
  ),
  filename        = NULL,
  category.names  = c("", ""),
  fill            = c("#E05A51", "#377EB8"),
  col             = c("#E05A51", "#377EB8"),
  alpha           = 0.4,
  cex             = 1.5,
  lwd             = 3,
  disable.logging = TRUE,
  fontfamily      = "sans",
  rotation.degree = -90,
  cat.pos         = c(0, 180),
  cat.dist        = c(0.04, 0.04)
)

# 3.3 ggplot and export
p <- ggplot() +
  annotation_custom(grobTree(venn.plot)) +
  theme_void() +
  theme(
    plot.title  = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.margin = margin(t = 10, r = 100, b = 10, l = 100)
  )

print(p)

ggsave("./Figure/Fenton_H2O2_venn_Fig3.svg",
       plot   = p,
       width  = 650/96,
       height = 400/96,
       device = "svg")