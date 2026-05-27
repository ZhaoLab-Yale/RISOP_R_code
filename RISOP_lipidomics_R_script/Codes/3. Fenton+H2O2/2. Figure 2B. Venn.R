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

group_info <- read_excel("./Input/Identification source.xlsx") %>%
  as.data.frame() %>%
  rename(Metabolite.curated = "Lipid", Group = "Separate.analysis") %>%
  select(Metabolite.curated, Group)

load(file = paste0("./Output/cmbd_lipids_Fenton.Rdata"))
load(file = paste0("./Output/cmbd_lipids_H2O2.Rdata"))

# ============================================================
# 2. Process lipid lists # -----
# ============================================================

cmbd_lipids_Fenton <- cmbd_lipids_Fenton %>%
  filter(grepl("<", Metabolite.curated, fixed = TRUE)) %>%
  select(Metabolite.curated) %>%
  mutate(Metabolite.curated = gsub("^(\\S+) (.+?>)\\s*(.*)", "\\1(\\2)\\3", Metabolite.curated))

cmbd_lipids_H2O2 <- cmbd_lipids_H2O2 %>%
  filter(grepl("<", Metabolite.curated, fixed = TRUE)) %>%
  select(Metabolite.curated) %>%
  mutate(Metabolite.curated = gsub("^(\\S+) (.+?>)\\s*(.*)", "\\1(\\2)\\3", Metabolite.curated))

union_metabolites <- union(cmbd_lipids_Fenton$Metabolite.curated, cmbd_lipids_H2O2$Metabolite.curated)

group_info <- group_info %>%
  filter(Metabolite.curated %in% union_metabolites)

# ============================================================
# 3. Create Venn diagram sets # -----
# ============================================================

# 3.1 Define sets by group
fenton_set <- group_info %>%
  filter(Group == "Fenton") %>%
  pull(Metabolite.curated)

h2o2_set <- group_info %>%
  filter(Group == "H2O2") %>%
  pull(Metabolite.curated)

both_set <- group_info %>%
  filter(Group == "Both") %>%
  pull(Metabolite.curated)

# 3.2 Assign "Both" lipids to both sets
fenton_all <- c(fenton_set, both_set)
h2o2_all   <- c(h2o2_set,   both_set)

# ============================================================
# 4. Plot and export Venn diagram # -----
# ============================================================

# 4.1 Draw Venn diagram
venn.plot <- venn.diagram(
  x = list(
    Fenton = fenton_all,
    H2O2   = h2o2_all
  ),
  filename         = NULL,
  category.names   = c("", ""),
  fill             = c("#E05A51", "#377EB8"),
  col              = c("#E05A51", "#377EB8"),
  alpha            = 0.4,
  cex              = 1.5,
  lwd              = 3,
  disable.logging  = TRUE,
  fontfamily       = "sans",
  rotation.degree  = -90,
  cat.pos          = c(0, 180),
  cat.dist         = c(0.04, 0.04)
)

# 4.2 Wrap in ggplot and export
p <- ggplot() +
  annotation_custom(grobTree(venn.plot)) +
  theme_void() +
  theme(plot.title  = element_text(hjust = 0.5, size = 18, face = "bold"),
        plot.margin = margin(t = 10, r = 100, b = 10, l = 100))

print(p)

ggsave("./Figure/Fenton_H2O2_venn.svg",
       plot   = p,
       width  = 650/96,
       height = 400/96,
       device = "svg")