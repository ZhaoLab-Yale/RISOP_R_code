rm(list = ls())

# ==============================================================================
# 1. Load required packages -----
# ==============================================================================
library(jsonlite)
library(dplyr)
library(stringr)
library(purrr)
library(readxl)

# ==============================================================================
# 2. Load input data -----
# ==============================================================================
intensity_mapping <- read_excel("./Input/ion intensity.xlsx")
data              <- fromJSON("./Input/KP4_All_PL.json")

# ==============================================================================
# 3. Helper functions -----
# ==============================================================================

# 3.1 Check if a fatty acid chain is oxidized (e.g., <OH>, <oxo>, <OOH>)
is_oxidized_chain <- function(fa_name) {
  return(str_detect(fa_name, "<[^>]*>"))
}

# 3.2 Look up ion intensity from mapping table by lipid name, adduct, and ion type
#     Name pattern examples: "PC" / "PC(O-" / "PC(P-"
#     Returns 999 if no match or intensity is NA
get_intensity <- function(lipid_name, adduct, ion_type, intensity_mapping) {
  if (str_detect(lipid_name, "\\(O-")) {
    name_pattern <- str_extract(lipid_name, "^[A-Z]+\\(O-")
  } else if (str_detect(lipid_name, "\\(P-")) {
    name_pattern <- str_extract(lipid_name, "^[A-Z]+\\(P-")
  } else {
    name_pattern <- str_extract(lipid_name, "^[A-Z]+")
  }
  
  match_row <- intensity_mapping %>%
    filter(name == name_pattern, Adduct == adduct)
  
  if (nrow(match_row) == 0) return(999)
  
  intensity <- switch(ion_type,
                      "precursor"             = match_row$`Precursor ion`[1],
                      "neutral_loss"          = match_row$`Neutral loss ion`[1],
                      "characteristic"        = match_row$`Characteristic ion`[1],
                      "oxidized_fragment"     = match_row$`Oxidized chain fragment`[1],
                      "non_oxidized_fragment" = match_row$`Non-oxidized chain fragment`[1],
                      999
  )
  
  if (is.na(intensity)) return(999)
  return(as.numeric(intensity))
}

# 3.3 Extract FA chain names from lipid name
#     Returns list(fa1, fa2); single-chain (lyso) lipids have fa2 = NULL
extract_fa_chains <- function(lipid_name) {
  if (str_detect(lipid_name, "\\([OP]-")) {
    # Ether lipids: e.g., PC(O-18:0_20:4) or PE(P-16:0_18:1)
    chains_match <- str_match(lipid_name, "\\(([OP]-\\S+?)_(\\S+?)\\)")
    if (!is.na(chains_match[1, 1])) {
      fa1 <- chains_match[1, 2]               # Keep O- or P- prefix as is
      fa2 <- paste0("FA", chains_match[1, 3])
      return(list(fa1 = fa1, fa2 = fa2))
    }
  } else {
    # Standard diacyl lipids
    chains_match <- str_match(lipid_name, "\\((\\S+?)_(\\S+?)\\)")
    if (!is.na(chains_match[1, 1])) {
      fa1 <- paste0("FA", chains_match[1, 2])
      fa2 <- paste0("FA", chains_match[1, 3])
      return(list(fa1 = fa1, fa2 = fa2))
    } else {
      # Lysophospholipid (single chain)
      chains_match <- str_match(lipid_name, "\\((\\S+?)\\)")
      if (!is.na(chains_match[1, 1])) {
        fa1 <- paste0("FA", chains_match[1, 2])
        return(list(fa1 = fa1, fa2 = NULL))
      }
    }
  }
  return(list(fa1 = NULL, fa2 = NULL))
}

# ==============================================================================
# 4. Main extraction function: parse one JSON entry into MSP fields -----
# ==============================================================================
extract_msp_data <- function(entry_name, entry_data, intensity_mapping) {
  
  sub_entry <- entry_data[[names(entry_data)[1]]]
  
  # 4.1 Basic metadata
  name           <- sub_entry[["name"]]
  precursor_mz   <- sub_entry[["mz"]]
  precursor_type <- sub_entry[["adduct"]]
  formula        <- sub_entry[["neutral_formula"]]
  
  origins <- sub_entry[["origins"]]
  lipid_origin <- if (!is.null(origins) && length(origins) > 0) {
    paste(origins, collapse = "; ")
  } else {
    ""
  }
  
  frag_info <- sub_entry[["ions"]][["frag_info"]]
  fa_chains <- extract_fa_chains(name)
  fa1_chain <- fa_chains$fa1
  fa2_chain <- fa_chains$fa2
  
  # 4.2 Derived metadata
  ion_mode       <- if (str_detect(precursor_type, "\\]-$")) "Negative" else "Positive"
  formatted_name <- str_replace(name, "([A-Z]+)\\((.+)\\)", "\\1 \\2")
  is_oxidized    <- str_detect(formatted_name, "<[^>]*>")
  is_ether_lipid <- str_detect(formatted_name, "[OP]-")
  
  # 4.3 Compound class (with Ether/Ox prefixes as needed)
  if (is_ether_lipid) {
    base_class  <- str_extract(formatted_name, "^[A-Z]+")
    lipid_class <- paste0("Ether", base_class)
  } else {
    lipid_class <- str_extract(formatted_name, "^[A-Z]+")
  }
  if (is_oxidized) lipid_class <- paste0("Ox", lipid_class)
  
  # 4.4 Electron-mass-corrected precursor m/z
  electron_mass          <- 0.00054858
  corrected_precursor_mz <- precursor_mz
  if (ion_mode == "Positive") {
    corrected_precursor_mz <- precursor_mz - electron_mass
  } else if (ion_mode == "Negative") {
    corrected_precursor_mz <- precursor_mz + electron_mass
  }
  
  # 4.5 Initialise ion tracking variables
  ion_list                      <- list()
  precursor_ion_mz              <- NULL
  neutral_loss_ion_mz           <- NULL
  characteristic_ion_mz         <- NULL
  oxidized_fragment_fa1_mzs     <- c()
  oxidized_fragment_fa2_mzs     <- c()
  non_oxidized_fragment_fa1_mzs <- c()
  non_oxidized_fragment_fa2_mzs <- c()
  all_fragment_fa1_mzs          <- c()
  all_fragment_fa2_mzs          <- c()
  
  # 4.6 Precursor ion
  precursor_intensity <- get_intensity(name, precursor_type, "precursor", intensity_mapping)
  precursor_ion_mz    <- round(corrected_precursor_mz, 4)
  ion_list[[as.character(precursor_ion_mz)]] <- list(
    mz        = corrected_precursor_mz,
    intensity = precursor_intensity,
    type      = "precursor"
  )
  
  # 4.7 Neutral loss ion (class- and adduct-specific)
  #     Check base class without "Ox" prefix
  base_lipid_class <- str_replace(lipid_class, "^Ox", "")
  special_nl <- NULL
  if ((base_lipid_class == "PC" || base_lipid_class == "LPC") && precursor_type == "[M+HCOO]-") {
    special_nl <- 60.0211294
  } else if (base_lipid_class == "EtherPC" && precursor_type == "[M+HCOO]-") {
    special_nl <- 60.0211294
  } else if (base_lipid_class == "PS" && precursor_type == "[M-H]-") {
    special_nl <- 87.032028
  }
  
  if (!is.null(special_nl)) {
    nl_mz        <- corrected_precursor_mz - special_nl
    nl_intensity <- get_intensity(name, precursor_type, "neutral_loss", intensity_mapping)
    neutral_loss_ion_mz <- round(nl_mz, 4)
    ion_list[[as.character(neutral_loss_ion_mz)]] <- list(
      mz        = nl_mz,
      intensity = nl_intensity,
      type      = "neutral_loss"
    )
  }
  
  # 4.8 Characteristic ion (head-group-specific)
  #     PC/LPC/PC O-: 224.0682 | PE/LPE/PE O-/PE P-: 196.0380
  #     PG/PS: 152.9958        | PI: 241.0119
  char_ion_mz <- NULL
  if      (str_detect(formatted_name, "^(PC|LPC|PC O-)"))       char_ion_mz <- 224.0682
  else if (str_detect(formatted_name, "^(PE|LPE|PE O-|PE P-)")) char_ion_mz <- 196.0380
  else if (str_detect(formatted_name, "^PG|PS"))                 char_ion_mz <- 152.9958
  else if (str_detect(formatted_name, "^PI"))                    char_ion_mz <- 241.0119
  
  if (!is.null(char_ion_mz)) {
    char_intensity        <- get_intensity(name, precursor_type, "characteristic", intensity_mapping)
    characteristic_ion_mz <- round(char_ion_mz, 4)
    ion_list[[as.character(characteristic_ion_mz)]] <- list(
      mz        = char_ion_mz,
      intensity = char_intensity,
      type      = "characteristic"
    )
  }
  
  # 4.9 Fragment ions from frag_info
  if (!is.null(frag_info) && length(frag_info) > 0) {
    for (frag_name in names(frag_info)) {
      frag    <- frag_info[[frag_name]]
      frag_mz <- frag[["mz"]]
      
      # Electron mass correction
      corrected_frag_mz <- frag_mz
      if (ion_mode == "Positive") {
        corrected_frag_mz <- frag_mz - electron_mass
      } else if (ion_mode == "Negative") {
        corrected_frag_mz <- frag_mz + electron_mass
      }
      
      chain_match <- str_match(frag_name, "^([^#]+)#")
      if (is.na(chain_match[1, 2])) next
      
      chain_name <- chain_match[1, 2]
      
      # Skip ether-linked (O-/P-) chains for ether lipids
      if (is_ether_lipid && str_detect(chain_name, "^[OP]-")) next
      
      # Assign fragment to FA1 or FA2
      is_fa1 <- !is.null(fa1_chain) && chain_name == fa1_chain
      is_fa2 <- !is.null(fa2_chain) && chain_name == fa2_chain
      
      # Intensity and oxidation tracking
      is_oxidized_frag <- is_oxidized_chain(chain_name)
      if (is_oxidized_frag) {
        intensity <- get_intensity(name, precursor_type, "oxidized_fragment", intensity_mapping)
        if (is_fa1) oxidized_fragment_fa1_mzs <- c(oxidized_fragment_fa1_mzs, round(corrected_frag_mz, 4))
        if (is_fa2) oxidized_fragment_fa2_mzs <- c(oxidized_fragment_fa2_mzs, round(corrected_frag_mz, 4))
      } else {
        intensity <- get_intensity(name, precursor_type, "non_oxidized_fragment", intensity_mapping)
        if (is_fa1) non_oxidized_fragment_fa1_mzs <- c(non_oxidized_fragment_fa1_mzs, round(corrected_frag_mz, 4))
        if (is_fa2) non_oxidized_fragment_fa2_mzs <- c(non_oxidized_fragment_fa2_mzs, round(corrected_frag_mz, 4))
      }
      
      if (is_fa1) all_fragment_fa1_mzs <- c(all_fragment_fa1_mzs, round(corrected_frag_mz, 4))
      if (is_fa2) all_fragment_fa2_mzs <- c(all_fragment_fa2_mzs, round(corrected_frag_mz, 4))
      
      # Add to ion list using m/z as key to avoid duplicates
      # Priority: precursor/NL/characteristic > fragment; for duplicate fragments keep lower intensity
      ion_key <- as.character(round(corrected_frag_mz, 4))
      if (!ion_key %in% names(ion_list)) {
        ion_list[[ion_key]] <- list(mz = corrected_frag_mz, intensity = intensity, type = "fragment")
      } else if (ion_list[[ion_key]][["type"]] == "fragment") {
        if (intensity < ion_list[[ion_key]][["intensity"]]) {
          ion_list[[ion_key]][["intensity"]] <- intensity
        }
      }
    }
  }
  
  return(list(
    NAME          = formatted_name,
    PRECURSORMZ   = corrected_precursor_mz,
    PRECURSORTYPE = precursor_type,
    FORMULA       = formula,
    IONMODE       = ion_mode,
    COMPOUNDCLASS = lipid_class,
    ION_LIST      = ion_list,
    ORIGIN        = lipid_origin,
    PRECURSOR_ION              = precursor_ion_mz,
    NEUTRAL_LOSS_ION           = neutral_loss_ion_mz,
    CHARACTERISTIC_ION         = characteristic_ion_mz,
    OXIDIZED_FRAGMENTS_FA1     = if (length(oxidized_fragment_fa1_mzs) > 0)     paste(oxidized_fragment_fa1_mzs,     collapse = ";") else NA,
    OXIDIZED_FRAGMENTS_FA2     = if (length(oxidized_fragment_fa2_mzs) > 0)     paste(oxidized_fragment_fa2_mzs,     collapse = ";") else NA,
    NON_OXIDIZED_FRAGMENTS_FA1 = if (length(non_oxidized_fragment_fa1_mzs) > 0) paste(non_oxidized_fragment_fa1_mzs, collapse = ";") else NA,
    NON_OXIDIZED_FRAGMENTS_FA2 = if (length(non_oxidized_fragment_fa2_mzs) > 0) paste(non_oxidized_fragment_fa2_mzs, collapse = ";") else NA,
    ALL_FRAGMENTS_FA1          = if (length(all_fragment_fa1_mzs) > 0)          paste(all_fragment_fa1_mzs,          collapse = ";") else NA,
    ALL_FRAGMENTS_FA2          = if (length(all_fragment_fa2_mzs) > 0)          paste(all_fragment_fa2_mzs,          collapse = ";") else NA
  ))
}

# ==============================================================================
# 5. Run extraction over all entries -----
# ==============================================================================
msp_data_list <- lapply(names(data), function(entry_name) {
  extract_msp_data(entry_name, data[[entry_name]], intensity_mapping)
})

# ==============================================================================
# 6. Build ion summary dataframe -----
# ==============================================================================
ion_df <- data.frame(
  Lipid_Name         = sapply(msp_data_list, function(x) x$NAME),
  Adduct             = sapply(msp_data_list, function(x) x$PRECURSORTYPE),
  Precursor_Ion      = sapply(msp_data_list, function(x) x$PRECURSOR_ION),
  Neutral_Loss_Ion   = sapply(msp_data_list, function(x) ifelse(is.null(x$NEUTRAL_LOSS_ION),    NA, x$NEUTRAL_LOSS_ION)),
  Characteristic_Ion = sapply(msp_data_list, function(x) ifelse(is.null(x$CHARACTERISTIC_ION), NA, x$CHARACTERISTIC_ION)),
  FA1                = sapply(msp_data_list, function(x) ifelse(is.na(x$ALL_FRAGMENTS_FA1),          NA, x$ALL_FRAGMENTS_FA1)),
  FA2                = sapply(msp_data_list, function(x) ifelse(is.na(x$ALL_FRAGMENTS_FA2),          NA, x$ALL_FRAGMENTS_FA2)),
  OxFA1              = sapply(msp_data_list, function(x) ifelse(is.na(x$OXIDIZED_FRAGMENTS_FA1),     NA, x$OXIDIZED_FRAGMENTS_FA1)),
  OxFA2              = sapply(msp_data_list, function(x) ifelse(is.na(x$OXIDIZED_FRAGMENTS_FA2),     NA, x$OXIDIZED_FRAGMENTS_FA2)),
  UnFA1              = sapply(msp_data_list, function(x) ifelse(is.na(x$NON_OXIDIZED_FRAGMENTS_FA1), NA, x$NON_OXIDIZED_FRAGMENTS_FA1)),
  UnFA2              = sapply(msp_data_list, function(x) ifelse(is.na(x$NON_OXIDIZED_FRAGMENTS_FA2), NA, x$NON_OXIDIZED_FRAGMENTS_FA2)),
  Origin             = sapply(msp_data_list, function(x) x$ORIGIN),
  stringsAsFactors = FALSE
)

# ==============================================================================
# 7. Convert to MSP format and write outputs -----
# ==============================================================================
msp_entries <- lapply(msp_data_list, function(entry) {
  ion_list <- entry[["ION_LIST"]]
  
  if (is.null(ion_list) || length(ion_list) == 0) {
    peaks     <- character(0)
    num_peaks <- 0
  } else {
    # Sort ions by m/z, then format as "mz\tintensity" lines
    ion_list_sorted <- ion_list[order(sapply(ion_list, function(x) x$mz))]
    peaks <- character(length(ion_list_sorted))
    for (i in seq_along(ion_list_sorted)) {
      ion      <- ion_list_sorted[[i]]
      peaks[i] <- paste(round(ion$mz, 4), ion$intensity, sep = "\t")
    }
    num_peaks <- length(ion_list_sorted)
  }
  
  entry_block <- c(
    paste("NAME:",          entry[["NAME"]]),
    paste("PRECURSORMZ:",   round(entry[["PRECURSORMZ"]], 4)),
    paste("PRECURSORTYPE:", entry[["PRECURSORTYPE"]]),
    "SMILES:",
    "INCHIKEY:",
    paste("FORMULA:",       entry[["FORMULA"]]),
    "RETENTIONTIME:",
    "CCS:",
    paste("IONMODE:",       entry[["IONMODE"]]),
    paste0("COMPOUNDCLASS: \"", entry[["COMPOUNDCLASS"]], "\""),
    paste("Comment:",       entry[["ORIGIN"]]),
    paste("Num Peaks:",     num_peaks),
    peaks
  )
  
  return(paste(entry_block, collapse = "\n"))
})

final_msp_content <- paste(msp_entries, collapse = "\n\n")
output_file       <- "./Output/KP4_All_PL.msp"

write.csv(ion_df, file = paste0("./Output/KP4_All_PL.csv"), row.names = FALSE)
writeLines(final_msp_content, output_file)

cat("Total entries:", length(msp_entries), "\n")