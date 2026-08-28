rm(list = ls())
library(tidyverse)
library(openxlsx)
library(stringr)
library(writexl)
library(readxl)
setwd(rstudioapi::getActiveProject())

# ==============================================================================
# 0. Shared functions -----
# ==============================================================================

# 0.1 Reformat lipid name: "PC 18:0_18:1" -> "PC(18:0_18:1)"
reformat_name <- function(x) {
  str_replace(x, "^(\\w+)\\s+(.+?)\\s*(\\([a-z]\\))?$", "\\1(\\2)\\3")
}

# 0.2 Reformat plasmalogen names to include both O- and P- forms
reformat_plasmalogen <- function(lipid_vec) {
  # capture chain part and trailing isomer suffix separately
  lpc_o <- str_match(lipid_vec, "^LPC\\(O-(\\d+):(\\d+)\\)(\\([a-z]\\))?$")
  lpe_p <- str_match(lipid_vec, "^LPE\\(P-(\\d+):(\\d+)\\)(\\([a-z]\\))?$")
  pe_o  <- str_match(lipid_vec, "^PE\\(O-(\\d+):(\\d+)_(.+?)\\)(\\([a-z]\\))?$")
  pc_o  <- str_match(lipid_vec, "^PC\\(O-(\\d+):(\\d+)_(.+?)\\)(\\([a-z]\\))?$")
  
  case_when(
    # ether chain must carry >=1 double bond to have a P- counterpart
    !is.na(lpc_o[,1]) & as.integer(lpc_o[,3]) > 0 ~ paste0(
      "LPC(O-", lpc_o[,2], ":", lpc_o[,3], ") or LPC(P-",
      lpc_o[,2], ":", as.integer(lpc_o[,3]) - 1, ")",
      replace_na(lpc_o[,4], "")
    ),
    !is.na(lpe_p[,1]) ~ paste0(
      "LPE(O-", lpe_p[,2], ":", as.integer(lpe_p[,3]) + 1, ") or LPE(P-",
      lpe_p[,2], ":", lpe_p[,3], ")",
      replace_na(lpe_p[,4], "")
    ),
    !is.na(pe_o[,1]) & as.integer(pe_o[,3]) > 0 ~ paste0(
      "PE(O-", pe_o[,2], ":", pe_o[,3], "_", pe_o[,4], ") or PE(P-",
      pe_o[,2], ":", as.integer(pe_o[,3]) - 1, "_", pe_o[,4], ")",
      replace_na(pe_o[,5], "")
    ),
    !is.na(pc_o[,1]) & as.integer(pc_o[,3]) > 0 ~ paste0(
      "PC(O-", pc_o[,2], ":", pc_o[,3], "_", pc_o[,4], ") or PC(P-",
      pc_o[,2], ":", as.integer(pc_o[,3]) - 1, "_", pc_o[,4], ")",
      replace_na(pc_o[,5], "")
    ),
    TRUE ~ lipid_vec
  )
}

# 0.3 Extract the standard matching name from a plasmalogen-reformatted lipid name
get_lipid_match <- function(lipid_vec) {
  part1 <- str_trim(str_split_fixed(lipid_vec, " or ", 2)[, 1])
  part2 <- str_trim(str_split_fixed(lipid_vec, " or ", 2)[, 2])
  suffix <- replace_na(str_extract(part2, "\\([a-z]\\)$"), "")
  case_when(
    # LPE original name is the P- form, which now sits in part2
    str_detect(part1, "^LPE\\(O-") ~ part2,
    # reattach the isomer suffix that now sits at the end of part2
    str_detect(part1, "^LPC\\(O-") ~ paste0(part1, suffix),
    str_detect(part1, "^PE\\(O-")  ~ paste0(part1, suffix),
    str_detect(part1, "^PC\\(O-")  ~ paste0(part1, suffix),
    TRUE                           ~ lipid_vec
  )
}

# 0.4 Imputation function: log2-transform, impute zeros from 5th-percentile normal, back-transform
impt.func <- function(norm.df) {
  norm.df <- norm.df[order(rownames(norm.df)), ]
  
  Norm_0_replace <- norm.df
  Norm_0_replace[Norm_0_replace == 0] <- NA
  log2df <- log2(Norm_0_replace)
  imptdf <- log2df
  
  set.seed(1234)
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

# ==============================================================================
# 1. Fenton -----
# ==============================================================================

# 1.1 Load data
df_curated     <- read.csv('./Input/df_curated_Fenton.csv', stringsAsFactors = F, skip = 0, header = T)
cust_curated   <- read.csv('./Input/cust_curated_Fenton.csv', stringsAsFactors = F, skip = 0, header = T)
Sample_info    <- read.csv('./Input/Sample_info.csv', stringsAsFactors = F, skip = 0, header = T)

# 1.2 Remove in-source fragmentation (ISF)
df_curated <- df_curated %>%
  filter(!Alignment.ID %in% c("38015", "39402", "30767", "39995", "41440", "39887"))

# 1.3 Normalize, flag, and remove rows with excessive zeros/NAs (threshold: 15/18 samples)
df_curated1 <- df_curated %>%
  select(-Polarity, -Adduct.type, -Average.Rt.min., -Ontology, -Alignment.ID) %>%
  mutate(across(where(is.numeric), ~ . / sum(., na.rm = TRUE))) %>%
  column_to_rownames(var = "Metabolite.curated") %>%
  mutate(Zero_percent = rowSums(across(everything(), ~ . == 0 | is.na(.))) / ncol(.))

zero_high_name_Fenton <- df_curated1 %>%
  filter(Zero_percent >= 15/18) %>%
  rownames() %>%
  reformat_name()

df_curated1 <- df_curated1 %>%
  filter(Zero_percent < 15/18) %>%
  select(-Zero_percent)

# 1.4 Apply imputation and back-transform from log2 scale
impt_df <- impt.func(df_curated1) %>%
  filter(rowSums(across(everything(), ~ . != 0)) > 0) %>%
  rownames_to_column(var = "Metabolite.curated") %>%
  mutate(across(where(is.numeric), ~ 2 ^ .))

# 1.5 Merge imputed data with Ontology annotation
cmbd_lipids_Fenton <- impt_df %>%
  left_join(df_curated %>% select(Metabolite.curated, Ontology), by = "Metabolite.curated")

# 1.6 Compute CV per timepoint
timepoints <- c("0h", "1min", "20min", "4h", "24h", "72h")

cv_df_Fenton <- timepoints %>%
  lapply(function(tp) {
    cols <- grep(paste0("Fenton_", tp, "_"), colnames(cmbd_lipids_Fenton), value = TRUE)
    cmbd_lipids_Fenton[, cols] %>%
      apply(1, function(x) {
        m <- mean(x)
        s <- sd(x)
        if (m == 0) return(NA)
        return(s / m * 100)
      }) %>%
      data.frame() %>%
      setNames(paste0("CV_", tp))
  }) %>%
  bind_cols(data.frame(Metabolite.curated = cmbd_lipids_Fenton$Metabolite.curated), .) %>%
  rowwise() %>%
  mutate(min_CV = if (!grepl("<", Metabolite.curated)) {
    min(c_across(c(CV_0h, CV_1min, CV_20min, CV_4h, CV_24h, CV_72h)), na.rm = TRUE)
  } else {
    min(c_across(c(CV_1min, CV_20min, CV_4h, CV_24h, CV_72h)), na.rm = TRUE)
  }) %>%
  ungroup()

# ==============================================================================
# 2. H2O2 -----
# ==============================================================================

# 2.1 Load data
df_curated     <- read.csv('./Input/df_curated_H2O2.csv', stringsAsFactors = F, skip = 0, header = T)
cust_curated   <- read.csv('./Input/cust_curated_H2O2.csv', stringsAsFactors = F, skip = 0, header = T)
Sample_info    <- read.csv('./Input/Sample_info.csv', stringsAsFactors = F, skip = 0, header = T)

# 2.2 Remove ISF
df_curated <- df_curated %>%
  filter(!Alignment.ID %in% c("38015", "39402", "30767", "39995", "41440", "39887"))

# 2.3 Normalize, flag, and remove rows with excessive zeros/NAs (threshold: 15/18 samples)
df_curated1 <- df_curated %>%
  select(-Polarity, -Adduct.type, -Average.Rt.min., -Ontology, -Alignment.ID) %>%
  mutate(across(where(is.numeric), ~ . / sum(., na.rm = TRUE))) %>%
  column_to_rownames(var = "Metabolite.curated") %>%
  mutate(Zero_percent = rowSums(across(everything(), ~ . == 0 | is.na(.))) / ncol(.))

zero_high_name_H2O2 <- df_curated1 %>%
  filter(Zero_percent >= 15/18) %>%
  rownames() %>%
  reformat_name()

df_curated1 <- df_curated1 %>%
  filter(Zero_percent < 15/18) %>%
  select(-Zero_percent)

# 2.4 Apply imputation and back-transform from log2 scale
impt_df <- impt.func(df_curated1) %>%
  filter(rowSums(across(everything(), ~ . != 0)) > 0) %>%
  rownames_to_column(var = "Metabolite.curated") %>%
  mutate(across(where(is.numeric), ~ 2 ^ .))

# 2.5 Merge imputed data with Ontology annotation
cmbd_lipids_H2O2 <- impt_df %>%
  left_join(df_curated %>% select(Metabolite.curated, Ontology), by = "Metabolite.curated")

# 2.6 Compute CV per timepoint
timepoints <- c("0h", "1min", "20min", "4h", "24h", "72h")

cv_df_H2O2 <- timepoints %>%
  lapply(function(tp) {
    cols <- grep(paste0("H2O2_", tp, "_"), colnames(cmbd_lipids_H2O2), value = TRUE)
    cmbd_lipids_H2O2[, cols] %>%
      apply(1, function(x) {
        m <- mean(x)
        s <- sd(x)
        if (m == 0) return(NA)
        return(s / m * 100)
      }) %>%
      data.frame() %>%
      setNames(paste0("CV_", tp))
  }) %>%
  bind_cols(data.frame(Metabolite.curated = cmbd_lipids_H2O2$Metabolite.curated), .) %>%
  rowwise() %>%
  mutate(min_CV = if (!grepl("<", Metabolite.curated)) {
    min(c_across(c(CV_0h, CV_1min, CV_20min, CV_4h, CV_24h, CV_72h)), na.rm = TRUE)
  } else {
    min(c_across(c(CV_1min, CV_20min, CV_4h, CV_24h, CV_72h)), na.rm = TRUE)
  }) %>%
  ungroup()

# ==============================================================================
# 3. Combined CV filtering -----
# ==============================================================================

# 3.1 Determine which lipids to remove
cv_both <- inner_join(
  cv_df_Fenton %>% select(Metabolite.curated, min_CV),
  cv_df_H2O2   %>% select(Metabolite.curated, min_CV),
  by = "Metabolite.curated",
  suffix = c("_Fenton", "_H2O2")
)

cv_Fenton_only <- cv_df_Fenton %>%
  filter(!Metabolite.curated %in% cv_df_H2O2$Metabolite.curated) %>%
  select(Metabolite.curated, min_CV)

cv_H2O2_only <- cv_df_H2O2 %>%
  filter(!Metabolite.curated %in% cv_df_Fenton$Metabolite.curated) %>%
  select(Metabolite.curated, min_CV)

cv_remove_Fenton <- c(
  cv_both        %>% filter(min_CV_Fenton >= 20 & min_CV_H2O2 >= 20) %>% pull(Metabolite.curated),
  cv_Fenton_only %>% filter(min_CV >= 20)                             %>% pull(Metabolite.curated)
)

cv_remove_H2O2 <- c(
  cv_both      %>% filter(min_CV_Fenton >= 20 & min_CV_H2O2 >= 20) %>% pull(Metabolite.curated),
  cv_H2O2_only %>% filter(min_CV >= 20)                             %>% pull(Metabolite.curated)
)

# 3.2 Apply combined CV filter
cmbd_lipids_Fenton <- cmbd_lipids_Fenton %>%
  filter(!Metabolite.curated %in% cv_remove_Fenton)

cmbd_lipids_H2O2 <- cmbd_lipids_H2O2 %>%
  filter(!Metabolite.curated %in% cv_remove_H2O2)

# 3.3 Save CV-removed names for supplementary table annotation
cv_removed_name_Fenton <- cv_remove_Fenton %>% reformat_name()
cv_removed_name_H2O2   <- cv_remove_H2O2   %>% reformat_name()

# 3.4 Save filtered datasets
save(cmbd_lipids_Fenton, file = paste0("./Output/cmbd_lipids_Fenton.Rdata"))
save(cmbd_lipids_H2O2,   file = paste0("./Output/cmbd_lipids_H2O2.Rdata"))

# ==============================================================================
# 4. Export supplementary tables -----
# ==============================================================================

# ------------------------------------------------------------------------------
# 4.1 Fenton -----
# ------------------------------------------------------------------------------

# Load raw Fenton data
df_curated <- read.csv('./Input/df_curated_Fenton.csv', stringsAsFactors = F, skip = 0, header = T)

# Reformat lipid names and rename columns for export
cmbd_lipids_for_supplementary <- df_curated %>%
  mutate(Metabolite.curated = reformat_name(Metabolite.curated)) %>%
  rename(Lipid = Metabolite.curated, `Retention time` = Average.Rt.min., Adduct = Adduct.type) %>%
  select(-Ontology, -Polarity) %>%
  relocate(Adduct, `Retention time`, .after = Lipid)

# Tab 1: Oxidized lipids, peak area (pre-filter, ISF not removed)
Ox_Fenton_prefiltered <- cmbd_lipids_for_supplementary %>%
  filter(grepl("<", Lipid, fixed = TRUE)) %>%
  mutate(Removal = if_else(Alignment.ID %in% c("38015", "39402", "30767", "39995", "41440", "39887"),
                           "ISF/misannotation", NA_character_)) %>%
  relocate(Removal, .after = Lipid)

# Tab 2: Non-oxidized lipids, peak area (pre-filter)
nonOx_Fenton_adduct_lookup <- cmbd_lipids_for_supplementary %>%
  filter(!grepl("<", Lipid, fixed = TRUE)) %>%
  select(Lipid, Adduct, `Retention time`)

nonOx_Fenton_prefiltered <- cmbd_lipids_for_supplementary %>%
  filter(!grepl("<", Lipid, fixed = TRUE)) %>%
  mutate(Removal = NA_character_, .after = Lipid) %>%
  mutate(Lipid = reformat_plasmalogen(Lipid))

# Tab 3: Oxidized lipids, normalized and filtered, with identification source
id_source <- read.xlsx("./Input/Identification source.xlsx", sheet = "Sheet1")

Ox_Fenton_filtered <- cmbd_lipids_Fenton %>%
  mutate(Metabolite.curated = reformat_name(Metabolite.curated)) %>%
  filter(grepl("<", Metabolite.curated, fixed = TRUE)) %>%
  rename("Lipid" = "Metabolite.curated") %>%
  left_join(Ox_Fenton_prefiltered %>% select(Lipid, Adduct, `Retention time`),
            by = "Lipid") %>%
  left_join(id_source %>% select(Lipid, Separate.analysis),
            by = "Lipid") %>%
  relocate(Adduct, `Retention time`, Separate.analysis, .after = Lipid) %>%
  select(-Ontology) %>%
  mutate(Separate.analysis = if_else(
    Separate.analysis == "Both",
    "Fenton, H2O2",
    Separate.analysis
  )) %>%
  rename("Identification in separate Fenton and H2O2-only references" = "Separate.analysis")

# Tab 4: Non-oxidized lipids, normalized and filtered
nonOx_Fenton_filtered <- cmbd_lipids_Fenton %>%
  mutate(Metabolite.curated = reformat_name(Metabolite.curated)) %>%
  filter(!grepl("<", Metabolite.curated, fixed = TRUE)) %>%
  rename("Lipid" = "Metabolite.curated") %>%
  left_join(nonOx_Fenton_adduct_lookup, by = "Lipid") %>%
  mutate(Lipid = reformat_plasmalogen(Lipid)) %>%
  relocate(Adduct, `Retention time`, .after = Lipid) %>%
  select(-Ontology)

# ------------------------------------------------------------------------------
# 4.2 H2O2 -----
# ------------------------------------------------------------------------------

# Load raw H2O2 data
df_curated <- read.csv('./Input/df_curated_H2O2.csv', stringsAsFactors = F, skip = 0, header = T)

# Reformat lipid names and rename columns for export
cmbd_lipids_for_supplementary <- df_curated %>%
  mutate(Metabolite.curated = reformat_name(Metabolite.curated)) %>%
  rename(Lipid = Metabolite.curated, `Retention time` = Average.Rt.min., Adduct = Adduct.type) %>%
  select(-Ontology, -Polarity) %>%
  relocate(Adduct, `Retention time`, .after = Lipid)

# Tab 1: Oxidized lipids, peak area (pre-filter, ISF not removed)
Ox_H2O2_prefiltered <- cmbd_lipids_for_supplementary %>%
  filter(grepl("<", Lipid, fixed = TRUE)) %>%
  mutate(Removal = if_else(Alignment.ID %in% c("38015", "39402", "30767", "39995", "41440", "39887"),
                           "ISF/misannotation", NA_character_)) %>%
  relocate(Removal, .after = Lipid)

# Tab 2: Non-oxidized lipids, peak area (pre-filter)
nonOx_H2O2_adduct_lookup <- cmbd_lipids_for_supplementary %>%
  filter(!grepl("<", Lipid, fixed = TRUE)) %>%
  select(Lipid, Adduct, `Retention time`)

nonOx_H2O2_prefiltered <- cmbd_lipids_for_supplementary %>%
  filter(!grepl("<", Lipid, fixed = TRUE)) %>%
  mutate(Removal = NA_character_, .after = Lipid) %>%
  mutate(Lipid = reformat_plasmalogen(Lipid))

# Tab 3: Oxidized lipids, normalized and filtered, with identification source
id_source <- read.xlsx("./Input/Identification source.xlsx", sheet = "Sheet1")

Ox_H2O2_filtered <- cmbd_lipids_H2O2 %>%
  mutate(Metabolite.curated = reformat_name(Metabolite.curated)) %>%
  filter(grepl("<", Metabolite.curated, fixed = TRUE)) %>%
  rename("Lipid" = "Metabolite.curated") %>%
  left_join(Ox_H2O2_prefiltered %>% select(Lipid, Adduct, `Retention time`),
            by = "Lipid") %>%
  left_join(id_source %>% select(Lipid, Separate.analysis),
            by = "Lipid") %>%
  relocate(Adduct, `Retention time`, Separate.analysis, .after = Lipid) %>%
  select(-Ontology) %>%
  mutate(Separate.analysis = if_else(
    Separate.analysis == "Both",
    "Fenton, H2O2",
    Separate.analysis
  )) %>%
  rename("Identification in separate Fenton and H2O2-only references" = "Separate.analysis")

# Tab 4: Non-oxidized lipids, normalized and filtered
nonOx_H2O2_filtered <- cmbd_lipids_H2O2 %>%
  mutate(Metabolite.curated = reformat_name(Metabolite.curated)) %>%
  filter(!grepl("<", Metabolite.curated, fixed = TRUE)) %>%
  rename("Lipid" = "Metabolite.curated") %>%
  left_join(nonOx_H2O2_adduct_lookup, by = "Lipid") %>%
  mutate(Lipid = reformat_plasmalogen(Lipid)) %>%
  relocate(Adduct, `Retention time`, .after = Lipid) %>%
  select(-Ontology)

# ------------------------------------------------------------------------------
# 4.3 Annotate removal reasons in pre-filtered tables -----
# ------------------------------------------------------------------------------

# Fenton - Ox
Ox_Fenton_prefiltered <- Ox_Fenton_prefiltered %>%
  mutate(
    Removal = if_else(Lipid %in% zero_high_name_Fenton  & is.na(Removal), "0 count too high", Removal),
    Removal = if_else(Lipid %in% cv_removed_name_Fenton & is.na(Removal), "CV>=20%",          Removal)
  ) %>%
  select(-Alignment.ID)

# Fenton - NonOx
nonOx_Fenton_prefiltered <- nonOx_Fenton_prefiltered %>%
  mutate(
    Lipid_match = get_lipid_match(Lipid),
    Removal = if_else(Lipid_match %in% zero_high_name_Fenton  & is.na(Removal), "0 count too high", Removal),
    Removal = if_else(Lipid_match %in% cv_removed_name_Fenton & is.na(Removal), "CV>=20%",          Removal)
  ) %>%
  select(-Lipid_match, -Alignment.ID)

# H2O2 - Ox
Ox_H2O2_prefiltered <- Ox_H2O2_prefiltered %>%
  mutate(
    Removal = if_else(Lipid %in% zero_high_name_H2O2  & is.na(Removal), "0 count too high", Removal),
    Removal = if_else(Lipid %in% cv_removed_name_H2O2 & is.na(Removal), "CV>=20%",          Removal)
  ) %>%
  select(-Alignment.ID)

# H2O2 - NonOx
nonOx_H2O2_prefiltered <- nonOx_H2O2_prefiltered %>%
  mutate(
    Lipid_match = get_lipid_match(Lipid),
    Removal = if_else(Lipid_match %in% zero_high_name_H2O2  & is.na(Removal), "0 count too high", Removal),
    Removal = if_else(Lipid_match %in% cv_removed_name_H2O2 & is.na(Removal), "CV>=20%",          Removal)
  ) %>%
  select(-Lipid_match, -Alignment.ID)

# ------------------------------------------------------------------------------
# 4.4 Add combined-analysis identification column -----
# ------------------------------------------------------------------------------
common_lipids <- intersect(Ox_Fenton_filtered$Lipid, Ox_H2O2_filtered$Lipid)

Ox_Fenton_filtered <- Ox_Fenton_filtered %>%
  mutate(`Identification in combined Fenton and H2O2-only references` = if_else(
    Lipid %in% common_lipids, "Fenton, H2O2", "Fenton"
  ), .after = `Identification in separate Fenton and H2O2-only references`)

Ox_H2O2_filtered <- Ox_H2O2_filtered %>%
  mutate(`Identification in combined Fenton and H2O2-only references` = if_else(
    Lipid %in% common_lipids, "Fenton, H2O2", "H2O2"
  ), .after = `Identification in separate Fenton and H2O2-only references`)

# ------------------------------------------------------------------------------
# 4.5 Export to Excel -----
# ------------------------------------------------------------------------------

list(
  `OxPLs (prefilter)`        = Ox_Fenton_prefiltered,
  `NonOx (prefilter)` = nonOx_Fenton_prefiltered,
  `OxPLs (filtered and normalized)`           = Ox_Fenton_filtered,
  `NonOx (filtered and normalized)`    = nonOx_Fenton_filtered
) %>%
  write.xlsx("./Output/Table S3. Intensities of OxPLs and non-oxidized lipids identified in Fenton references.xlsx")

list(
  `OxPLs (prefilter)`        = Ox_H2O2_prefiltered,
  `NonOx (prefilter)` = nonOx_H2O2_prefiltered,
  `OxPLs (filtered and normalized)`           = Ox_H2O2_filtered,
  `NonOx (filtered and normalized)`    = nonOx_H2O2_filtered
) %>%
  write.xlsx("./Output/Table S4. Intensities of OxPLs and non-oxidized lipids identified in H2O2 references.xlsx")

