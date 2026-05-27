rm(list = ls())
library(tidyverse)
library(openxlsx)
library(stringr)
library(readxl)
library(writexl)
setwd(rstudioapi::getActiveProject())

# ============================================================
# 1. Load data # -----
# ============================================================

Norm_method = "Protein_norm"

df_curated     <- read.csv('./Input/df_curated_ML210.csv',    stringsAsFactors = F, skip = 0, header = T)
cust_curated   <- read.csv('./Input/cust_curated_ML210.csv',  stringsAsFactors = F, skip = 0, header = T)
species_for_QT <- read.csv('./Input/species_for_QT.csv',      stringsAsFactors = F, skip = 0, header = T)
Sample_info    <- read.csv('./Input/Sample_info.csv',         stringsAsFactors = F, skip = 0, header = T)

# ============================================================
# 2. Pre-process input data # -----
# ============================================================

df_curated <- df_curated %>%
  filter(!Alignment.ID %in% c("38015", "39402", "30767", "39995", "41440", "39887", "42220")) %>%
  select(-Alignment.ID)

df_curated <- df_curated %>%
  mutate(New_ontology = str_remove_all(Ontology, "Ox|Ether"))

cust_curated <- cust_curated %>%
  mutate(New_ontology = str_remove_all(Ontology, "Ox|Ether"))

sample_vector <- Sample_info %>%
  filter(is.na(Outlier) | Outlier != "Yes") %>%
  pull(Sample_name) %>%
  intersect(colnames(df_curated))

species_for_QT <- species_for_QT %>%
  mutate(
    Stock.cont.uM      = as.numeric(Stock.cont.uM),
    Volume.used.uL     = as.numeric(Volume.used.uL),
    Volume.for.recs.uL = as.numeric(Volume.for.recs.uL),
    IS.Cont            = Stock.cont.uM * Volume.used.uL / Volume.for.recs.uL
  )

# ============================================================
# 3. Reshape to long format # -----
# ============================================================

df_long <- df_curated %>%
  select(Average.Rt.min.,
         Metabolite.curated,
         Adduct.type,
         Ontology,
         New_ontology,
         sample_vector) %>%
  pivot_longer(
    cols      = sample_vector,
    names_to  = "Sample_name",
    values_to = "Area"
  )

IS_long <- cust_curated %>%
  filter(Ontology %in% species_for_QT$Ontology & Adduct.type %in% species_for_QT$Adduct) %>%
  select(Adduct.type, Ontology, sample_vector) %>%
  left_join(species_for_QT, by = c("Adduct.type", "Ontology")) %>%
  pivot_longer(
    cols      = sample_vector,
    names_to  = "Sample_name",
    values_to = "IS_area"
  )

all_lipid_info <- df_long %>%
  left_join(IS_long %>%
              select(Adduct.type, Ontology, Sample_name, IS.Cont, IS_area),
            by = c("Ontology", "Sample_name", "Adduct.type")) %>%
  mutate(IS_area = ifelse(!is.na(IS_area), IS_area, NA))

# ============================================================
# 4. Quantification # -----
# ============================================================

lipid_abs <- all_lipid_info %>%
  mutate(abs_cont = Area / IS_area * IS.Cont) %>%
  na.omit()

Lipid_rela <- all_lipid_info[rowSums(is.na(all_lipid_info)) > 0, ] %>%
  .[, colSums(is.na(.)) == 0]

# ============================================================
# 5. Compute normalization factors # -----
# ============================================================

all_mean_median <- lipid_abs %>%
  group_by(Sample_name) %>%
  summarize(
    total_mean   = mean(abs_cont,   na.rm = TRUE),
    total_median = median(abs_cont, na.rm = TRUE),
    total_sum    = sum(abs_cont,    na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    NF_total_mean   = total_mean   / mean(total_mean),
    NF_total_median = total_median / median(total_median),
    NF_total_sum    = total_sum    / mean(total_sum)
  )

PL_mean_median <- lipid_abs %>%
  filter(Ontology %in% c("PA", "PC", "PE", "PI", "PS", "PG", "LPC", "LPE")) %>%
  group_by(Sample_name) %>%
  summarize(
    PL_mean   = mean(abs_cont,   na.rm = TRUE),
    PL_median = median(abs_cont, na.rm = TRUE),
    PL_sum    = sum(abs_cont,    na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    NF_PL_mean   = PL_mean   / mean(PL_mean),
    NF_PL_median = PL_median / median(PL_median),
    NF_PL_sum    = PL_sum    / mean(PL_sum)
  )

Class_area_sum <- all_lipid_info %>%
  group_by(New_ontology, Adduct.type, Sample_name) %>%
  summarize(
    Class_peak_sum = sum(Area, na.rm = TRUE),
    .groups = "drop"
  )

Total_area <- all_lipid_info %>%
  group_by(Sample_name) %>%
  summarize(
    Total_peak_sum = sum(Area, na.rm = TRUE),
    .groups = "drop"
  )

Class_cont_sum <- lipid_abs %>%
  group_by(New_ontology, Adduct.type, Sample_name) %>%
  summarize(
    Class_concern_sum = sum(abs_cont, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 6. Merge normalization factors and apply normalization # -----
# ============================================================

val_for_abs_norm <- lipid_abs %>%
  left_join(all_mean_median, by = "Sample_name") %>%
  left_join(PL_mean_median,  by = "Sample_name") %>%
  left_join(Class_area_sum,  by = c("Sample_name", "New_ontology", "Adduct.type")) %>%
  left_join(Total_area,      by = "Sample_name") %>%
  left_join(Class_cont_sum,  by = c("Sample_name", "New_ontology", "Adduct.type")) %>%
  left_join(Sample_info,     by = "Sample_name")

val_for_rela_norm <- Lipid_rela %>%
  left_join(all_mean_median, by = "Sample_name") %>%
  left_join(PL_mean_median,  by = "Sample_name") %>%
  left_join(Class_area_sum,  by = c("Sample_name", "New_ontology", "Adduct.type")) %>%
  left_join(Total_area,      by = "Sample_name") %>%
  left_join(Class_cont_sum,  by = c("Sample_name", "New_ontology", "Adduct.type")) %>%
  left_join(Sample_info,     by = "Sample_name")

lip_abs_norm <- val_for_abs_norm %>%
  mutate(
    all_median_norm   = abs_cont / NF_total_median,
    all_mean_norm     = abs_cont / NF_total_mean,
    all_sum_norm      = abs_cont / NF_total_sum,
    PL_median_norm    = abs_cont / NF_PL_median,
    PL_mean_norm      = abs_cont / NF_PL_mean,
    PL_sum_norm       = abs_cont / NF_PL_sum,
    Protein_norm      = (abs_cont / Protein) * 150,
    Percent_area_norm = Area / Class_peak_sum,
    Percent_cont_norm = abs_cont / Class_concern_sum,
    Total_area_norm   = Area / Total_peak_sum,
    No                = abs_cont
  )

lip_rela_norm <- val_for_rela_norm %>%
  mutate(
    all_median_norm   = Area / NF_total_median,
    all_mean_norm     = Area / NF_total_mean,
    all_sum_norm      = Area / NF_total_sum,
    PL_median_norm    = Area / NF_PL_median,
    PL_mean_norm      = Area / NF_PL_mean,
    PL_sum_norm       = Area / NF_PL_sum,
    Protein_norm      = Area / Protein,
    Percent_area_norm = Area / Class_peak_sum,
    Percent_cont_norm = Area / Class_peak_sum,
    Total_area_norm   = Area / Total_peak_sum,
    No                = Area
  )

# ============================================================
# 7. Pivot wide, filter low-coverage rows, and impute # -----
# ============================================================

abs_norm_wd <- lip_abs_norm %>%
  select(c("Metabolite.curated", "Sample_name", !!sym(Norm_method))) %>%
  pivot_wider(
    names_from  = Sample_name,
    values_from = !!sym(Norm_method)
  ) %>%
  column_to_rownames(var = "Metabolite.curated") %>%
  mutate(Zero_percent = rowSums(across(everything(), ~ . == 0 | is.na(.))) / ncol(.))

rela_norm_wd <- lip_rela_norm %>%
  select(c("Metabolite.curated", "Sample_name", !!sym(Norm_method))) %>%
  pivot_wider(
    names_from  = Sample_name,
    values_from = !!sym(Norm_method)
  ) %>%
  column_to_rownames(var = "Metabolite.curated") %>%
  mutate(Zero_percent = rowSums(across(everything(), ~ . == 0 | is.na(.))) / ncol(.))

abs_norm_wd <- abs_norm_wd %>%
  filter(Zero_percent < 15/18) %>%
  select(-Zero_percent)

rela_norm_wd <- rela_norm_wd %>%
  filter(Zero_percent < 15/18) %>%
  select(-Zero_percent)

impt.func <- function(norm.df) {
  Norm_0_replace <- norm.df
  Norm_0_replace[Norm_0_replace == 0] <- NA
  log2df <- log2(Norm_0_replace)
  imptdf <- log2df
  
  set.seed(12345)
  for (i in 1:nrow(imptdf)) {
    b <- imptdf[i, which(is.na(imptdf[i, ]) == T)]
    
    imputed_values <- rnorm(n    = length(b),
                            mean = quantile(as.numeric(imptdf[i, ]), na.rm = T, 0.05),
                            sd   = abs(0.3 * quantile(as.numeric(imptdf[i, ]), na.rm = T, 0.05)))
    
    max_allowed    <- min(imptdf[i, ], na.rm = T)
    imputed_values <- pmin(imputed_values, max_allowed)
    
    imptdf[i, which(is.na(imptdf[i, ]) == T)] <- imputed_values
  }
  return(imptdf)
}

impt_abs <- impt.func(abs_norm_wd) %>%
  filter(rowSums(across(everything(), ~ . != 0)) > 0) %>%
  rownames_to_column(var = "Metabolite.curated") %>%
  left_join(df_long %>% select(c("Metabolite.curated", "Ontology")) %>% distinct(),
            by = "Metabolite.curated") %>%
  mutate(across(where(is.numeric), ~ 2 ^ .))

impt_rela <- impt.func(rela_norm_wd) %>%
  filter(rowSums(across(everything(), ~ . != 0)) > 0) %>%
  rownames_to_column(var = "Metabolite.curated") %>%
  left_join(df_long %>% select(c("Metabolite.curated", "Ontology")) %>% distinct(),
            by = "Metabolite.curated") %>%
  mutate(across(where(is.numeric), ~ 2 ^ .))

# ============================================================
# 8. Combined lipid table # -----
# ============================================================

cmbd_lipids <- rbind(impt_abs)

cmbd_lipids <- cmbd_lipids %>%
  mutate(Ontology = case_when(
    str_detect(Metabolite.curated, "-") ~ {
      letter_before_dash <- str_extract(Metabolite.curated, "[A-Za-z](?=-)")
      paste0(Ontology, "-", letter_before_dash) %>%
        str_remove_all("Ether")
    },
    TRUE ~ Ontology
  ))

# ============================================================
# 9. CV-based filtering # -----
# ============================================================

timepoints <- c("0min", "20min", "45min", "80min", "120min", "270min")

cv_df <- timepoints %>%
  lapply(function(tp) {
    cols <- grep(paste0("ML210_", tp, "_"), colnames(cmbd_lipids), value = TRUE)
    cmbd_lipids[, cols] %>%
      apply(1, function(x) {
        m <- mean(x)
        s <- sd(x)
        if (m == 0) return(NA)
        return(s / m * 100)
      }) %>%
      data.frame() %>%
      setNames(paste0("CV_", tp))
  }) %>%
  bind_cols(data.frame(Metabolite.curated = cmbd_lipids$Metabolite.curated), .) %>%
  rowwise() %>%
  mutate(min_CV = if (!grepl("<", Metabolite.curated)) {
    min(c_across(c(CV_0min, CV_20min, CV_45min, CV_80min, CV_120min, CV_270min)), na.rm = TRUE)
  } else {
    min(c_across(c(CV_20min, CV_45min, CV_80min, CV_120min, CV_270min)), na.rm = TRUE)
  }) %>%
  ungroup()

cmbd_lipids <- cmbd_lipids %>%
  filter(Metabolite.curated %in% (cv_df %>% filter(min_CV < 20) %>% pull(Metabolite.curated)))

# ============================================================
# 10. Save filtered lipid table # -----
# ============================================================

save(cmbd_lipids, file = paste0("./Output/cmbd_lipids.Rdata"))

# ============================================================
# 11. Build supplementary tables and export # -----
# ============================================================

format_lipid_names <- function(df) {
  df %>%
    mutate(Metabolite.curated = gsub("^(\\w+) ([^ (]+)( \\(([a-z])\\))?$",
                                     "\\1(\\2)\\3",
                                     Metabolite.curated) %>%
             gsub(" \\(([a-z])\\)", "(\\1)", .)) %>%
    rename(Lipid = Metabolite.curated,
           `Retention time` = Average.Rt.min.,
           Adduct = Adduct.type) %>%
    select(-Ontology) %>%
    relocate(Adduct, `Retention time`, .after = Lipid)
}

rename_ether_lipids <- function(df) {
  df %>%
    mutate(Lipid = case_when(
      str_detect(Lipid, "^LPC\\(O-\\d+:\\d+\\)$") ~ paste0(
        Lipid, " or LPC(P-",
        str_match(Lipid, "LPC\\(O-(\\d+):(\\d+)\\)")[,2], ":",
        as.integer(str_match(Lipid, "LPC\\(O-(\\d+):(\\d+)\\)")[,3]) - 1, ")"
      ),
      str_detect(Lipid, "^LPE\\(P-\\d+:\\d+\\)$") ~ paste0(
        "LPE(O-",
        str_match(Lipid, "LPE\\(P-(\\d+):(\\d+)\\)")[,2], ":",
        as.integer(str_match(Lipid, "LPE\\(P-(\\d+):(\\d+)\\)")[,3]) + 1,
        ") or ", Lipid
      ),
      TRUE ~ Lipid
    ))
}

# 11.3 Prefiltered supplementary table (before CV filter)
cmbd_lipids_pre <- rbind(impt_abs) %>%
  mutate(Ontology = case_when(
    str_detect(Metabolite.curated, "-") ~ {
      letter_before_dash <- str_extract(Metabolite.curated, "[A-Za-z](?=-)")
      paste0(Ontology, "-", letter_before_dash) %>%
        str_remove_all("Ether")
    },
    TRUE ~ Ontology
  )) %>%
  left_join(df_curated %>% select(Metabolite.curated, Adduct.type, Average.Rt.min.),
            by = "Metabolite.curated") %>%
  format_lipid_names()

cv_df_formatted <- cv_df %>%
  mutate(Metabolite.curated = gsub("^(\\w+) ([^ (]+)( \\(([a-z])\\))?$",
                                   "\\1(\\2)\\3",
                                   Metabolite.curated) %>%
           gsub(" \\(([a-z])\\)", "(\\1)", .))

# OxPLs prefiltered：
Ox_ML210_prefiltered <- cmbd_lipids_pre %>%
  filter(grepl("<", Lipid, fixed = TRUE)) %>%
  left_join(cv_df_formatted %>% select(Metabolite.curated, min_CV),
            by = c("Lipid" = "Metabolite.curated")) %>%
  mutate(Removal = ifelse(min_CV >= 20, "CV>=20%", NA)) %>%
  select(-min_CV) %>%
  relocate(Removal, .after = Lipid)

# NonOx prefiltered：
NonOx_ML210_prefiltered <- cmbd_lipids_pre %>%
  filter(!grepl("<", Lipid, fixed = TRUE)) %>%
  left_join(cv_df_formatted %>% select(Metabolite.curated, min_CV),
            by = c("Lipid" = "Metabolite.curated")) %>%
  mutate(Removal = ifelse(min_CV >= 20, "CV>=20%", NA)) %>%
  select(-min_CV) %>%
  relocate(Removal, .after = Lipid) %>%
  rename_ether_lipids()

# 11.4 Filtered supplementary table (after CV filter)
id_source <- read_xlsx('./Input/Identification source.xlsx', sheet = 2)

cmbd_lipids_post <- cmbd_lipids %>%
  left_join(df_curated %>% select(Metabolite.curated, Adduct.type, Average.Rt.min.),
            by = "Metabolite.curated") %>%
  format_lipid_names()

Ox_ML210_filtered <- cmbd_lipids_post %>%
  filter(grepl("<", Lipid, fixed = TRUE)) %>%
  left_join(id_source %>% select(Metabolite.curated, Group),
            by = c("Lipid" = "Metabolite.curated")) %>%
  relocate(Adduct, `Retention time`, Group, .after = Lipid) %>%
  rename(`Identification source` = Group)

NonOx_ML210_filtered <- cmbd_lipids_post %>%
  filter(!grepl("<", Lipid, fixed = TRUE)) %>%
  rename_ether_lipids()

# 11.5 Export
list(
  `OxPLs (prefiltered)`        = Ox_ML210_prefiltered,
  `NonOx_lipids (prefiltered)` = NonOx_ML210_prefiltered,
  `OxPLs (filtered)`           = Ox_ML210_filtered,
  `NonOx_lipids (filtered)`    = NonOx_ML210_filtered
) %>%
  write.xlsx("./Output/Table S5_ML210.xlsx")
