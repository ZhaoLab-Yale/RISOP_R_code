rm(list = ls())
library(tidyverse)
library(ggplot2)
library(ggtext)
setwd(rstudioapi::getActiveProject())

# ==============================================================================
# 1. Load input data -----
# ==============================================================================
cust_curated <- read.csv('./Input/cust_curated_ML210_QC.csv', stringsAsFactors = F, skip = 0, header = T)
Sample_info  <- read.csv('./Input/Sample_info.csv',           stringsAsFactors = F, skip = 0, header = T)

# ==============================================================================
# 2. Data filtering and preparation -----
# ==============================================================================

# 2.1 Keep only ML210_QC samples
samples_to_keep <- Sample_info %>%
  filter(Group == "ML210_QC") %>%
  pull(Sample_name)

# 2.2 Filter by valid ontology/adduct combinations and reformat metabolite names
cust_curated <- cust_curated %>%
  select(1:5, all_of(samples_to_keep)) %>%
  filter(
    (Ontology == "PC"  & Adduct.type == "[M+H]+")   |
      (Ontology == "PE"  & Adduct.type == "[M-H]-")   |
      (Ontology == "LPE" & Adduct.type == "[M-H]-")   |
      (Ontology == "CL"  & Adduct.type == "[M-H]-")   |
      (Ontology == "FA"  & Adduct.type == "[M-H]-")   |
      (Ontology == "DG"  & Adduct.type == "[M+NH4]+") |
      (Ontology == "TG"  & Adduct.type == "[M+NH4]+") |
      (Ontology == "PI"  & Adduct.type == "[M-H]-")   |
      (Ontology == "PS"  & Adduct.type == "[M-H]-")   |
      (Ontology == "PG"  & Adduct.type == "[M-H]-")   |
      (Ontology == "LPC" & Adduct.type == "[M+H]+")
  ) %>%
  mutate(Metabolite.curated = gsub("^(\\w+) (.+)$", "\\1(\\2)", Metabolite.curated))

# 2.3 Reshape to long format
df_long <- cust_curated %>%
  pivot_longer(
    cols      = all_of(samples_to_keep),
    names_to  = "Sample",
    values_to = "Intensity"
  )

# 2.4 Normalize each lipid by its group mean (mean = 1 per lipid)
df_normalized <- df_long %>%
  group_by(Metabolite.curated) %>%
  mutate(Normalized_Intensity = Intensity / mean(Intensity, na.rm = TRUE)) %>%
  ungroup()

# ==============================================================================
# 3. Summary statistics -----
# ==============================================================================

# 3.1 Helper function: format a number as HTML scientific notation
format_sci_html <- function(x) {
  exp  <- floor(log10(abs(x)))
  coef <- x / 10^exp
  sprintf("%.2f \u00d7 10<sup>%d</sup>", coef, exp)
}

# 3.2 Use the first lipid (CL) to display "CV:" and "Int:" prefix labels
first_lipid <- unique(cust_curated$Metabolite.curated)[1]

# 3.3 Calculate mean, SEM, CV, and raw mean per lipid
#     The first lipid entry receives "CV:" and "Int:" prefixes in the annotation label
df_summary <- df_normalized %>%
  group_by(Metabolite.curated) %>%
  summarise(
    mean_int = mean(Normalized_Intensity, na.rm = TRUE),
    sem_int  = sd(Normalized_Intensity,   na.rm = TRUE) / sqrt(n()),
    cv_pct   = (sd(Normalized_Intensity,  na.rm = TRUE) / mean(Normalized_Intensity, na.rm = TRUE)) * 100,
    mean_raw = mean(Intensity, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(
    Metabolite.curated = factor(Metabolite.curated,
                                levels = unique(cust_curated$Metabolite.curated)),
    label_html = ifelse(
      Metabolite.curated == first_lipid,
      paste0(
        "<span style='color:#3651AD'>CV: ", round(cv_pct, 1), "%</span><br>",
        "<span style='color:#784814'>Int: ", sapply(mean_raw, format_sci_html), "</span>"
      ),
      paste0(
        "<span style='color:#3651AD'>", round(cv_pct, 1), "%</span><br>",
        "<span style='color:#784814'>", sapply(mean_raw, format_sci_html), "</span>"
      )
    )
  )

# ==============================================================================
# 4. QC plot -----
# ==============================================================================

p <- ggplot(df_summary, aes(x = Metabolite.curated, y = mean_int)) +
  # Error bars (SEM)
  geom_errorbar(aes(ymin = mean_int - sem_int,
                    ymax = mean_int + sem_int),
                width = 0.3, colour = "black") +
  # Mean crossbar
  geom_crossbar(aes(ymin = mean_int, ymax = mean_int),
                width = 0.4, colour = "#1F4E79", linewidth = 0.5) +
  # Individual replicate points
  geom_jitter(data = df_normalized,
              aes(x = factor(Metabolite.curated, levels = levels(df_summary$Metabolite.curated)),
                  y = Normalized_Intensity),
              width = 0.15, size = 3, alpha = 0.9, colour = "#4C78A8") +
  # CV and raw intensity annotations above each bar
  geom_richtext(aes(y = mean_int + sem_int + 0.03,
                    label = label_html),
                vjust = -0.8, size = 4.2,
                colour = "black", lineheight = 1.4,
                fill = NA, label.color = NA) +
  # ±10% CV reference lines
  geom_hline(yintercept = 1.1, linetype = "dashed", colour = "#787878", linewidth = 0.4) +
  geom_hline(yintercept = 0.9, linetype = "dashed", colour = "#787878", linewidth = 0.4) +
  annotate("text", x = Inf, y = 1.1, label = "+10% CV", hjust = 1.1, vjust = -0.8, size = 4.2, colour = "#787878") +
  annotate("text", x = Inf, y = 0.9, label = "-10% CV", hjust = 1.1, vjust =  1.8, size = 4.2, colour = "#787878") +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.12)),
    limits = c(0.88, 1.10),
    breaks = scales::pretty_breaks()
  ) +
  scale_x_discrete(
    expand = expansion(add = c(0.7, 0.7))
  ) +
  labs(
    x = "Internal standard",
    y = "Relative intensity"
  ) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 13),
    axis.text.y      = element_text(size = 13),
    axis.title.x     = element_text(size = 18, face = "bold", margin = margin(t = -12, b = 0)),
    axis.title.y     = element_text(size = 18, face = "bold", margin = margin(r = 16)),
    plot.margin      = margin(t = 10, r = 10, b = 0, l = 40, unit = "pt"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p)

# ==============================================================================
# 5. Export -----
# ==============================================================================

pdf("./Figure/QC_lipid_intensity.pdf", width = 1100/96, height = 700/96)
print(p)
dev.off()