rm(list = ls())
library(tidyverse)
library(ggrepel)
library(readxl)
library(dplyr)
library(stringr)
library(ggbreak)
library(gridExtra)
library(grid)
library(cowplot)
library(scales)
library(purrr)
setwd(rstudioapi::getActiveProject())

# ==============================================================================
# 1. Load input data -----
# ==============================================================================

NonOx_df      <- read.csv('./Input/NonOx_df.csv',      stringsAsFactors = F, skip = 0, header = T)
df_curated_Ox <- read.csv('./Input/df_curated_Ox.csv', stringsAsFactors = F, skip = 0, header = T)

Ox_df_curated <- df_curated_Ox %>%
  arrange(Metabolite.curated, Average.Rt.min.) %>%
  group_by(Metabolite.curated) |>
  mutate(
    Metabolite.curated = if (n() > 1) {
      paste0(Metabolite.curated, " (", letters[1:n()][row_number()], ")")
    } else {
      Metabolite.curated
    }
  ) |>
  ungroup()

ISF_IDs <- c("38015", "39402", "30767", "39995", "41440", "39887")
Ox_df_curated <- Ox_df_curated %>%
  mutate(is_ISF = Alignment.ID %in% ISF_IDs)

df_lib <- read.csv('./Output/KP4_All_PL.csv', stringsAsFactors = F, skip = 0, header = T)

# ==============================================================================
# 2. Combine non-oxidized and oxidized lipids -----
# ==============================================================================
NonOx_df_flagged <- NonOx_df %>%
  mutate(is_ISF = FALSE)

df <- NonOx_df_flagged %>%
  bind_rows(Ox_df_curated) %>%
  select(intersect(names(NonOx_df_flagged), names(Ox_df_curated)))

# ==============================================================================
# 3. Separate oxidized and non-oxidized subsets -----
# ==============================================================================
df_Ox <- df %>%
  filter(grepl("<", Metabolite.curated, fixed = TRUE)) %>%
  select(Metabolite.curated, Average.Rt.min., Average.Mz, is_ISF) %>%
  mutate(Metabolite_wo_abc = str_replace(Metabolite.curated, "\\s*\\([a-z]\\)$", "")) %>%
  left_join(df_lib[, c("Lipid_Name", "Origin")], by = c("Metabolite_wo_abc" = "Lipid_Name"))

df_nonOx <- df %>%
  filter(!grepl("<", Metabolite.curated, fixed = TRUE)) %>%
  select(Metabolite.curated, Average.Rt.min., Average.Mz, is_ISF)

# ==============================================================================
# 4. Build combined lipid table with precursor linkage -----
# ==============================================================================
df_all <- df_Ox %>%
  filter(!is.na(Origin))

original_NonOx <- df_all %>%
  pull(Origin) %>%
  strsplit(";") %>%
  unlist() %>%
  trimws() %>%
  unique() %>%
  gsub("\\(", " ", .) %>%
  gsub("\\)", "", .)

original_NonOx_true <- df_nonOx %>%
  filter(Metabolite.curated %in% original_NonOx)

all_lipids <- df_all %>%
  bind_rows(original_NonOx_true) %>%
  arrange(Metabolite.curated, Average.Rt.min.) %>%
  group_by(Metabolite.curated) |>
  mutate(
    Metabolite.curated = if (n() > 1) {
      paste0(Metabolite.curated, " (", letters[1:n()][row_number()], ")")
    } else {
      Metabolite.curated
    }
  ) |>
  ungroup()

all_lipids <- all_lipids %>%
  mutate(Metabolite.curated = str_replace(Metabolite.curated,
                                          "^(\\w+)\\s+(.+?)\\s*(\\([a-z]\\))?$", "\\1(\\2)\\3")) %>%
  mutate(Group = Origin) %>%
  mutate(Group = if_else(
    is.na(Group) | Group == "",
    str_extract(Metabolite.curated, "^[^(]+\\([^)]+\\)"),
    Group
  ))

all_lipids <- all_lipids %>%
  separate_rows(Group, sep = "; ")

all_lipids <- all_lipids %>%
  group_by(Group) %>%
  filter(any(str_detect(Metabolite.curated, fixed(Group)))) %>%
  ungroup()

# ==============================================================================
# 5. Identify truncated lipids -----
# ==============================================================================
all_lipids2 <- all_lipids %>%
  mutate(Group2 = Group) %>%
  mutate(Group2 = ifelse(
    str_detect(Metabolite.curated, "_") &
      sapply(Metabolite.curated, function(x) {
        matches <- str_extract_all(x, "(\\d+):")[[1]]
        if (length(matches) >= 2) {
          num1 <- as.numeric(str_extract(matches[1], "\\d+"))
          num2 <- as.numeric(str_extract(matches[2], "\\d+"))
          return(num1 > num2)
        }
        return(FALSE)
      }),
    paste0(Group2, "_truncated"),
    Group2
  ))

all_lipids2 <- all_lipids2 %>%
  mutate(Group3 = Group2) %>%
  mutate(Group3 = ifelse(
    grepl("_truncated", Group3),
    gsub("<[a-zA-Z]+>", "", Metabolite_wo_abc) %>%
      gsub("([A-Za-z]+) ([0-9])", "\\1(\\2", .) %>%
      paste0(")"),
    Group3
  ))

all_lipids3 <- all_lipids2 %>%
  {
    truncated_rows <- filter(., grepl("_truncated", Group2))
    
    copied_rows <- truncated_rows %>%
      rowwise() %>%
      do({
        current_group  <- .$Group
        current_group3 <- .$Group3
        
        all_lipids2 %>%
          filter(grepl(current_group, Metabolite.curated, fixed = TRUE)) %>%
          mutate(Group3 = current_group3)
      }) %>%
      ungroup()
    
    bind_rows(., copied_rows)
  } %>%
  distinct()

# ==============================================================================
# 6. Kendrick Mass Defect (KMD) settings -----
# ==============================================================================
base_nominal <- 16.000000
base_exact   <- 15.994914620
scale_factor <- base_nominal / base_exact

KMD_CLUSTER_GAP <- 0.005

# ==============================================================================
# 7. Color and shape mappings -----
# ==============================================================================
my_fill_colors <- c("Non-oxidized"    = "#EE30A7",
                    "Oxidized_Group1" = "#AB82FF",
                    "Oxidized_Group2" = "#FFB5C5",
                    "ISF"             = "darkgreen")

my_outline_colors <- c("Non-oxidized"    = "#A5006E",
                       "Oxidized_Group1" = "#6A3FBF",
                       "Oxidized_Group2" = "#D97090",
                       "ISF"             = "grey30")

shape_values <- c(21, 22, 23, 24)
ISF_SHAPE    <- 25

# ==============================================================================
# 8. Legend grob builder -----
# ==============================================================================
make_legend_grob <- function(panel_data,
                             dot_size     = 4,
                             text_size    = 7,
                             title_size   = text_size,
                             row_height   = 0.4,
                             units        = "npc",
                             x_dot        = 0.03,
                             gap_dot_num  = 0.08,
                             gap_num_name = 0.06,
                             dot_alpha    = 0.8,
                             isf_note     = "(potential ISF product, removed)",
                             title_ox     = "Annotated OxPLs",
                             title_isf    = "Misannotated/ISF",
                             title_nonox  = "Precursor lipid",
                             entry_gap    = 0,
                             title_gap    = 0,
                             section_gap  = 1,
                             y_npc        = 0) {
  
  carbondb_values <- unique(panel_data$CarbonDB_for_shape[panel_data$color_group != "ISF"])
  carbondb_values <- carbondb_values[!is.na(carbondb_values)]
  
  carbondb_sorted <- c(
    intersect("matched", carbondb_values),
    setdiff(carbondb_values, "matched") %>% sort()
  )
  
  pch_lookup <- setNames(shape_values[seq_along(carbondb_sorted)], carbondb_sorted)
  
  legend_df <- panel_data %>%
    select(point_label, color_group, Metabolite.curated, CarbonDB_for_shape) %>%
    distinct()
  
  ox_df <- legend_df %>%
    filter(!color_group %in% c("Non-oxidized", "ISF")) %>%
    mutate(modification = str_extract(Metabolite.curated, "(?<=<)[^>]+(?=>)")) %>%
    group_by(modification) %>%
    mutate(mod_order = min(point_label)) %>%
    ungroup() %>%
    arrange(mod_order, point_label) %>%
    select(-modification, -mod_order)
  
  isf_df <- legend_df %>%
    filter(color_group == "ISF")
  
  nonox_df <- legend_df %>%
    filter(color_group == "Non-oxidized")
  
  sections <- list(
    list(title = title_ox,    df = ox_df),
    list(title = title_isf,   df = isf_df),
    list(title = title_nonox, df = nonox_df)
  )
  sections <- Filter(function(s) nrow(s$df) > 0, sections)
  
  items <- list()
  for (k in seq_along(sections)) {
    s <- sections[[k]]
    if (k > 1 && section_gap > 0) {
      items[[length(items) + 1]] <- list(kind = "spacer", weight = section_gap)
    }
    items[[length(items) + 1]] <- list(kind = "title", label = s$title, weight = 1)
    if (title_gap > 0) {
      items[[length(items) + 1]] <- list(kind = "spacer", weight = title_gap)
    }
    for (i in seq_len(nrow(s$df))) {
      row_i <- s$df[i, ]
      if (i > 1 && entry_gap > 0) {
        items[[length(items) + 1]] <- list(kind = "spacer", weight = entry_gap)
      }
      items[[length(items) + 1]] <- list(
        kind   = "entry",
        row    = row_i,
        weight = 1
      )
    }
  }
  
  total_rows      <- sum(sapply(items, function(it) it$weight))
  total_height_cm <- total_rows * row_height
  
  x_num  <- x_dot + gap_dot_num
  x_name <- x_num + gap_num_name
  ux     <- function(v) unit(v, units)
  
  grob_list        <- list()
  cumulative_units <- 0
  
  for (it in items) {
    w <- it$weight
    
    row_center_units <- cumulative_units + w / 2
    y_pos            <- 1 - row_center_units / total_rows
    cumulative_units <- cumulative_units + w
    
    if (it$kind == "spacer") next
    
    if (it$kind == "title") {
      grob_list[[length(grob_list) + 1]] <- textGrob(
        label = it$label,
        x = ux(x_dot) - unit(dot_size / 2, "pt"),
        y = unit(y_pos, "npc"),
        just = "left",
        gp = gpar(fontsize = title_size, fontface = "bold", col = "black"))
      next
    }
    
    row_i        <- it$row
    is_isf_row   <- row_i$color_group == "ISF"
    is_nonox_row <- row_i$color_group == "Non-oxidized"
    
    if (is_isf_row) {
      pch_val <- ISF_SHAPE
    } else {
      pch_val <- if (!is.na(row_i$CarbonDB_for_shape) &&
                     row_i$CarbonDB_for_shape %in% names(pch_lookup)) {
        pch_lookup[[row_i$CarbonDB_for_shape]]
      } else {
        21
      }
    }
    
    fill_col    <- alpha(my_fill_colors[row_i$color_group], dot_alpha)
    outline_col <- my_outline_colors[row_i$color_group]
    
    if (is_isf_row || is_nonox_row) {
      grob_list[[length(grob_list) + 1]] <- pointsGrob(
        x = ux(x_dot), y = unit(y_pos, "npc"),
        pch = pch_val, size = unit(dot_size, "pt"),
        gp = gpar(fill = fill_col, col = outline_col, lwd = 0.8))
      grob_list[[length(grob_list) + 1]] <- textGrob(
        label = row_i$Metabolite.curated,
        x = ux(x_num), y = unit(y_pos, "npc"),
        just = "left", gp = gpar(fontsize = text_size, col = "black"))
    } else {
      grob_list[[length(grob_list) + 1]] <- pointsGrob(
        x = ux(x_dot), y = unit(y_pos, "npc"),
        pch = pch_val, size = unit(dot_size, "pt"),
        gp = gpar(fill = fill_col, col = outline_col, lwd = 0.8))
      grob_list[[length(grob_list) + 1]] <- textGrob(
        label = as.character(row_i$point_label),
        x = ux(x_num), y = unit(y_pos, "npc"),
        just = "left", gp = gpar(fontsize = text_size, col = "black"))
      grob_list[[length(grob_list) + 1]] <- textGrob(
        label = row_i$Metabolite.curated,
        x = ux(x_name), y = unit(y_pos, "npc"),
        just = "left", gp = gpar(fontsize = text_size, col = "black"))
    }
  }
  
  gTree(
    children = do.call(gList, grob_list),
    vp = viewport(
      width  = unit(1, "npc"),
      height = unit(total_height_cm, "cm"),
      y      = unit(y_npc, "npc"),
      just   = c("center", "bottom")
    )
  )
}

# ==============================================================================
# 8b. Measure the width a legend actually needs
# ==============================================================================
legend_width_cm <- function(panel_data,
                            text_size,
                            title_size,
                            x_dot        = 0.10,
                            gap_dot_num  = 0.30,
                            gap_num_name = 0.32,
                            title_nonox  = "Precursor lipid",
                            pad_cm       = 0.15) {
  
  x_num  <- x_dot + gap_dot_num
  x_name <- x_num + gap_num_name
  
  tw <- function(lbl, fs, bold = FALSE) {
    if (length(lbl) == 0) return(0)
    max(vapply(lbl, function(l) {
      convertWidth(
        grobWidth(textGrob(l, gp = gpar(
          fontsize = fs,
          fontface = if (bold) "bold" else "plain"))),
        "cm", valueOnly = TRUE)
    }, numeric(1)))
  }
  
  legend_df <- panel_data %>%
    select(point_label, color_group, Metabolite.curated) %>%
    distinct()
  
  ox_names    <- legend_df$Metabolite.curated[!legend_df$color_group %in% c("Non-oxidized", "ISF")]
  plain_names <- legend_df$Metabolite.curated[ legend_df$color_group %in% c("Non-oxidized", "ISF")]
  
  w_ox    <- if (length(ox_names))    x_name + tw(ox_names,    text_size) else 0
  w_plain <- if (length(plain_names)) x_num  + tw(plain_names, text_size) else 0
  
  titles   <- c("Annotated OxPLs", "Misannotated/ISF", title_nonox)
  w_titles <- x_dot + tw(titles, title_size, bold = TRUE)
  
  max(w_ox, w_plain, w_titles) + pad_cm
}

# ==============================================================================
# 8c. Strip white rect backgrounds from a grob tree (for ggbreak panels) -----
# ==============================================================================
make_grob_transparent <- function(g) {
  if (inherits(g, "rect")) {
    g$gp$fill <- "transparent"
    g$gp$col  <- NA
    return(g)
  }
  if (!is.null(g$children)) {
    g$children <- do.call(gList, lapply(g$children, make_grob_transparent))
  }
  if (!is.null(g$grobs)) {
    g$grobs <- lapply(g$grobs, make_grob_transparent)
  }
  g
}

# ==============================================================================
# 9. Compute KMD values and assign color groups -----
# ==============================================================================
For_mass_defect <- all_lipids3 %>%
  mutate(
    Kendrick_mass         = Average.Mz * scale_factor,
    Nominal_Kendrick_mass = round(Kendrick_mass),
    KMD                   = Kendrick_mass - Nominal_Kendrick_mass,
    is_oxidized           = grepl("<", Metabolite.curated),
    precursor             = gsub("<.*$", "", Metabolite.curated) %>%
      gsub("\\s*\\([a-z]\\)\\s*$", "", .) %>%
      trimws(),
    oxidation_status      = factor(is_oxidized,
                                   levels = c(FALSE, TRUE),
                                   labels = c("Non-oxidized", "Oxidized"))
  )

For_mass_defect <- For_mass_defect %>%
  group_by(Group3) %>%
  mutate(KMD_rounded = {
    ord     <- order(KMD)
    cl      <- integer(n())
    cl[ord] <- cumsum(c(1, as.integer(diff(KMD[ord]) > KMD_CLUSTER_GAP)))
    cl
  }) %>%
  ungroup() %>%
  group_by(Group3) %>%
  mutate(
    KMD_row_id = ifelse(oxidation_status == "Oxidized",
                        as.numeric(factor(KMD_rounded,
                                          levels = sort(unique(KMD_rounded[oxidation_status == "Oxidized"])))),
                        NA)
  ) %>%
  ungroup()

For_mass_defect <- For_mass_defect %>%
  group_by(Group3) %>%
  mutate(n_unique_kmd = length(unique(KMD_row_id[!is.na(KMD_row_id)]))) %>%
  ungroup() %>%
  group_by(Group3, KMD_row_id) %>%
  mutate(
    min_RT_in_row = ifelse(oxidation_status == "Oxidized", min(Average.Rt.min.), NA)
  ) %>%
  ungroup()

For_mass_defect <- For_mass_defect %>%
  group_by(Group3) %>%
  mutate(
    KMD_row_id_reordered = ifelse(
      n_unique_kmd >= 2 & oxidation_status == "Oxidized",
      as.numeric(factor(min_RT_in_row, levels = sort(unique(min_RT_in_row[!is.na(min_RT_in_row)])))),
      KMD_row_id
    )
  ) %>%
  ungroup()

For_mass_defect <- For_mass_defect %>%
  group_by(Group3) %>%
  mutate(
    color_group = case_when(
      is_ISF                                           ~ "ISF",
      oxidation_status == "Non-oxidized"               ~ "Non-oxidized",
      n_unique_kmd == 1                                ~ "Oxidized_Group2",
      n_unique_kmd >= 2 & (KMD_row_id_reordered == 1) ~ "Oxidized_Group2",
      n_unique_kmd >= 2 & (KMD_row_id_reordered == 2) ~ "Oxidized_Group1",
      TRUE                                             ~ "Other"
    )
  ) %>%
  ungroup()

For_mass_defect <- For_mass_defect %>%
  distinct(Group3, Metabolite.curated, Average.Rt.min., Average.Mz, .keep_all = TRUE)

For_mass_defect <- For_mass_defect %>%
  group_by(Group3) %>%
  mutate(
    kmd_group_key = ifelse(
      oxidation_status == "Non-oxidized",
      as.character(KMD_rounded),
      paste0("ox_", KMD_row_id_reordered)
    )
  ) %>%
  ungroup() %>%
  group_by(Group3, kmd_group_key) %>%
  mutate(min_RT_per_kmd_group = min(Average.Rt.min.)) %>%
  ungroup() %>%
  group_by(Group3) %>%
  mutate(
    kmd_group_order = as.numeric(factor(
      min_RT_per_kmd_group,
      levels = sort(unique(min_RT_per_kmd_group))
    ))
  ) %>%
  arrange(Group3, kmd_group_order, Average.Rt.min.) %>%
  mutate(point_label = row_number()) %>%
  select(-kmd_group_key, -min_RT_per_kmd_group, -kmd_group_order,
         -min_RT_in_row, -KMD_row_id_reordered) %>%
  ungroup()

For_mass_defect <- For_mass_defect %>%
  group_by(Group3) %>%
  mutate(point_label = ifelse(color_group == "ISF",
                              NA_integer_,
                              cumsum(color_group != "ISF"))) %>%
  ungroup()

For_mass_defect <- For_mass_defect %>%
  group_by(Group3) %>%
  filter(any(oxidation_status == "Oxidized")) %>%
  ungroup()

For_mass_defect <- For_mass_defect %>%
  mutate(
    CarbonDB = case_when(
      oxidation_status == "Oxidized"     ~ str_extract(Metabolite.curated, "(?<=\\()[0-9:_]+(?=<)"),
      oxidation_status == "Non-oxidized" ~ str_extract(Metabolite.curated, "(?<=\\()[0-9:_]+(?=\\))"),
      TRUE                               ~ NA_character_
    )
  )

For_mass_defect <- For_mass_defect %>%
  mutate(Group3_CarbonDB = str_extract(Group3, "(?<=\\()[0-9:_]+(?=\\))")) %>%
  group_by(Group3) %>%
  mutate(
    CarbonDB_for_shape = case_when(
      color_group == "ISF"        ~ NA_character_,
      CarbonDB == Group3_CarbonDB ~ "matched",
      TRUE                        ~ CarbonDB
    )
  ) %>%
  select(-Group3_CarbonDB) %>%
  ungroup()

# ==============================================================================
# 10. Label positioning function -----
# ==============================================================================
position_labels_by_kmd <- function(sub_df,
                                   relative_label_width  = 0.05,
                                   relative_spacing      = 0.08,
                                   relative_y_offset     = 0.06,
                                   cluster_gap_threshold = 0.15) {
  
  sub_df <- as.data.frame(sub_df)
  sub_df$.row_id <- seq_len(nrow(sub_df))
  rt_col  <- sub_df[["Average.Rt.min."]]
  kmd_col <- sub_df[["KMD"]]
  
  x_range <- max(rt_col,  na.rm = TRUE) - min(rt_col,  na.rm = TRUE)
  y_range <- max(kmd_col, na.rm = TRUE) - min(kmd_col, na.rm = TRUE)
  
  absolute_label_width <- x_range * relative_label_width
  absolute_spacing     <- x_range * relative_spacing
  absolute_y_offset    <- y_range * relative_y_offset
  absolute_cluster_gap <- x_range * cluster_gap_threshold
  
  sub_df$KMD_group <- round(kmd_col / 0.02) * 0.02
  
  label_positions <- sub_df %>%
    as.data.frame() %>%
    group_by(KMD_group) %>%
    arrange(across(all_of("Average.Rt.min."))) %>%
    mutate(
      rt_diff     = c(0, diff(.data[["Average.Rt.min."]])),
      sub_cluster = cumsum(rt_diff > absolute_cluster_gap) + 1
    ) %>%
    group_by(KMD_group, sub_cluster) %>%
    mutate(
      label_y         = first(.data[["KMD"]]) + absolute_y_offset,
      n_in_subcluster = n()
    ) %>%
    ungroup()
  
  label_positions <- label_positions %>%
    group_by(KMD_group, sub_cluster) %>%
    mutate(
      label_x = {
        n       <- first(n_in_subcluster)
        rt_vals <- .data[["Average.Rt.min."]]
        if (n == 1) {
          rt_vals
        } else {
          rt_sorted <- sort(rt_vals)
          min_dist  <- min(diff(rt_sorted))
          if (min_dist >= absolute_label_width) {
            rt_vals
          } else {
            center_x         <- mean(rt_vals)
            required_spacing <- (n - 1) * absolute_spacing
            half_span        <- required_spacing / 2
            seq(center_x - half_span, center_x + half_span, length.out = n)
          }
        }
      }
    ) %>%
    ungroup() %>%
    select(.row_id, label_x, label_y)
  
  sub_df %>%
    left_join(label_positions, by = ".row_id") %>%
    select(-.row_id)
}

# ==============================================================================
# 11. Plot tuning parameters -----
# ==============================================================================

# 11a. GRID set
GRID_DOT_SIZE    <- 6
GRID_ROW_HEIGHT  <- 0.35
GRID_ENTRY_GAP   <- 0
GRID_TITLE_GAP   <- 0
GRID_X_DOT       <- 0.10
GRID_GAP_DN      <- 0.30
GRID_GAP_NN      <- 0.32
GRID_DOT_ALPHA   <- 0.8
GRID_TICK_LENGTH <- 0.05
GRID_SECTION_GAP <- 0.4
GRID_PT_SIZE     <- 2

GRID_PANEL_W_CM <- 4.0
GRID_LEG_PAD_CM <- -1

GRID_TEXT_SIZE       <- 7.5
GRID_LEG_TITLE_SIZE  <- 7.5
GRID_TITLE_SIZE      <- 9
GRID_AXIS_TITLE_SIZE <- 9
GRID_AXIS_TEXT_SIZE  <- 8
GRID_PT_LABEL_SIZE   <- 8

GRID_N_XBREAKS  <- 3
GRID_XBREAK_PAD <- 0.10
GRID_LEG_Y_NPC  <- -0
GRID_MARGIN_L   <- 2

# 11b. SOLO set
SOLO_DOT_SIZE    <- 6
SOLO_ROW_HEIGHT  <- 0.35
SOLO_ENTRY_GAP   <- 0
SOLO_TITLE_GAP   <- 0
SOLO_GAP_DN      <- 0.025
SOLO_GAP_NN      <- 0.04
SOLO_DOT_ALPHA   <- 0.8
SOLO_LEG_WIDTH   <- 1.1
SOLO_TICK_LENGTH <- 0.05
SOLO_SECTION_GAP <- 0.4
SOLO_PT_SIZE     <- 2

SOLO_TEXT_SIZE       <- 8.3
SOLO_LEG_TITLE_SIZE  <- 9
SOLO_TITLE_SIZE      <- 9.5
SOLO_AXIS_TITLE_SIZE <- 10
SOLO_AXIS_TEXT_SIZE  <- 9
SOLO_PT_LABEL_SIZE   <- 8

SOLO_LEG_Y_NPC     <- -0
SOLO_YTITLE_MARGIN <- -5

# 11c. Shared annotation settings
ISF_NOTE_MAIN  <- "(misannotation, removed)"
ISF_NOTE_BREAK <- "(potential ISF product, removed)"

ISF_LINE_WIDTH_MIN <- 0
ISF_LINE_Y_OFFSET  <- 0.003
ISF_LINE_TYPE      <- "12"
ISF_LINE_WIDTH     <- 0.3

# 11d. Per-group x-axis range / break settings
xlim_settings <- list(
  "PC(16:0_16:1)" = c(9.4, 14),
  "PC(16:0_18:1)" = c(9.9, 21),
  "PC(16:0_22:5)" = c(12, 13.5),
  "PS(18:0_20:3)" = c(12.95, 13.23),
  "PS(18:1_20:1)" = c(11, 17),
  "PS(18:1_22:1)" = c(15, 16)
)

xbreak_settings <- list(
  "PC(16:0_16:1)" = list(gap = c(10.7, 11), scale_ratio = 0.4),
  "PC(16:0_18:1)" = list(gap = c(11.3, 13), scale_ratio = 0.4)
)

# ==============================================================================
# 12. Panel-building function (GRID set) -----
# ==============================================================================
For_mass_defect <- For_mass_defect %>%
  group_by(Group3) %>%
  filter(any(color_group == "ISF")) %>%
  ungroup()

build_panel <- function(grp, break_range = NULL, scale_ratio = "fixed") {
  
  sub_df <- For_mass_defect %>%
    filter(Group3 == grp) %>%
    as.data.frame()
  
  sub_df <- position_labels_by_kmd(sub_df,
                                   relative_label_width  = 0.05,
                                   relative_spacing      = 0.08,
                                   relative_y_offset     = 0.06,
                                   cluster_gap_threshold = 0.15)
  
  sub_df <- sub_df %>%
    mutate(RT_cluster = as.integer(factor(round(Average.Rt.min. / 0.05) * 0.05))) %>%
    group_by(RT_cluster) %>%
    mutate(label_y = max(KMD) + 0.0003) %>%
    ungroup() %>%
    select(-RT_cluster)
  
  sub_df <- sub_df %>%
    mutate(precursor_label = ifelse(color_group == "Non-oxidized",
                                    str_extract(Metabolite.curated, "\\([a-z]\\)$"),
                                    NA_character_))
  
  x_rng <- if (grp %in% names(xlim_settings)) {
    xlim_settings[[grp]]
  } else {
    range(sub_df$Average.Rt.min., na.rm = TRUE)
  }
  x_inset        <- diff(x_rng) * GRID_XBREAK_PAD
  x_breaks_panel <- seq(x_rng[1] + x_inset, x_rng[2] - x_inset, length.out = GRID_N_XBREAKS)
  
  carbondb_values  <- unique(sub_df$CarbonDB_for_shape[sub_df$color_group != "ISF"])
  carbondb_values  <- carbondb_values[!is.na(carbondb_values)]
  carbondb_sorted  <- c(intersect("matched", carbondb_values),
                        setdiff(carbondb_values, "matched") %>% sort())
  shape_vals_panel <- setNames(shape_values[seq_along(carbondb_sorted)], carbondb_sorted)
  
  sub_df <- sub_df %>%
    mutate(CarbonDB_for_shape = ifelse(color_group == "ISF", "ISF_marker", CarbonDB_for_shape))
  shape_vals_panel["ISF_marker"] <- ISF_SHAPE
  
  colors_used          <- unique(sub_df$color_group)
  fill_colors_panel    <- my_fill_colors[colors_used]
  outline_colors_panel <- my_outline_colors[colors_used]
  
  n_precursors   <- sum(sub_df$oxidation_status == "Non-oxidized")
  precursor_word <- ifelse(n_precursors > 1, "precursors", "precursor")
  panel_title    <- paste0("Ox", grp, " and ", precursor_word)
  
  p_base <- ggplot(sub_df,
                   aes(x = Average.Rt.min., y = KMD,
                       fill = color_group, color = color_group, shape = CarbonDB_for_shape)) +
    geom_segment(
      data = subset(sub_df, !(color_group == "Non-oxidized" & is.na(precursor_label)) & color_group != "ISF"),
      aes(x = Average.Rt.min., y = KMD, xend = label_x, yend = label_y),
      color = "grey40", size = 0.3) +
    geom_point(size = GRID_PT_SIZE, alpha = GRID_DOT_ALPHA, stroke = 0.5) +
    geom_text(data = subset(sub_df, color_group != "Non-oxidized" & !is.na(point_label)),
              aes(x = label_x, y = label_y, label = point_label),
              color = "black", size = GRID_PT_LABEL_SIZE / .pt, hjust = 0.5, vjust = -0.2) +
    geom_text(data = subset(sub_df, color_group == "Non-oxidized" & !is.na(precursor_label)),
              aes(x = label_x, y = label_y, label = precursor_label),
              color = "black", size = GRID_PT_LABEL_SIZE / .pt, hjust = 0.5, vjust = -0.2) +
    scale_fill_manual(values = fill_colors_panel) +
    scale_color_manual(values = outline_colors_panel) +
    scale_shape_manual(values = shape_vals_panel) +
    scale_x_continuous(breaks = x_breaks_panel,
                       labels = scales::label_number(accuracy = 0.01)) +
    scale_y_continuous(expand = expansion(mult = 0.25)) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title          = element_text(face = "bold", size = GRID_TITLE_SIZE, color = "black",
                                         hjust = 0.5, margin = margin(b = 2, unit = "pt")),
      axis.title.x        = element_text(face = "bold", size = GRID_AXIS_TITLE_SIZE, color = "black"),
      axis.title.y        = element_text(face = "bold", size = GRID_AXIS_TITLE_SIZE, color = "black"),
      axis.text.x         = element_text(size = GRID_AXIS_TEXT_SIZE, color = "black",
                                         margin = margin(t = 8, unit = "pt")),
      axis.text.y         = element_text(size = GRID_AXIS_TEXT_SIZE, color = "black",
                                         margin = margin(r = 8, unit = "pt")),
      legend.position     = "none",
      panel.grid.major    = element_blank(),
      panel.grid.minor    = element_blank(),
      panel.border        = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.ticks          = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length.x = unit(-GRID_TICK_LENGTH, "cm"),
      axis.ticks.length.y = unit(-GRID_TICK_LENGTH, "cm"),
      plot.margin         = margin(t = 4, r = 4, b = 4, l = GRID_MARGIN_L),
      plot.background     = element_rect(fill = "transparent", color = NA),
      panel.background    = element_rect(fill = "transparent", color = NA)
    )
  
  p <- p_base + labs(title = panel_title, x = "Retention time (min)", y = "O-KMD")
  
  if (is.null(break_range)) {
    if (grp %in% names(xlim_settings)) {
      p <- p + coord_cartesian(xlim = xlim_settings[[grp]])
    }
  } else {
    p <- p +
      theme(panel.border = element_blank(),
            axis.line     = element_line(color = "black", linewidth = 0.5))
    if (grp %in% names(xlim_settings)) {
      p <- p + coord_cartesian(xlim = xlim_settings[[grp]])
    }
    p <- p + scale_x_break(break_range, scales = scale_ratio)
    p <- ggdraw() + draw_grob(grid.grabExpr(print(p)))
  }
  
  sub_df_leg <- sub_df %>%
    mutate(CarbonDB_for_shape = ifelse(color_group == "ISF", NA_character_, CarbonDB_for_shape)) %>%
    mutate(
      Metabolite.curated = ifelse(color_group == "Non-oxidized",
                                  str_remove(Metabolite.curated, "\\([a-z]\\)$"),
                                  Metabolite.curated),
      point_label = ifelse(color_group == "Non-oxidized", NA_integer_, point_label)
    ) %>%
    group_by(Metabolite.curated) %>%
    mutate(
      Metabolite_leg = ifelse(
        color_group == "Non-oxidized" & any(!is.na(precursor_label)),
        paste0(Metabolite.curated,
               paste(sort(unique(precursor_label[!is.na(precursor_label)])), collapse = "/")),
        Metabolite.curated
      )
    ) %>%
    ungroup() %>%
    mutate(Metabolite.curated = Metabolite_leg) %>%
    select(-Metabolite_leg)
  
  leg <- make_legend_grob(sub_df_leg,
                          dot_size     = GRID_DOT_SIZE,
                          text_size    = GRID_TEXT_SIZE,
                          title_size   = GRID_LEG_TITLE_SIZE,
                          row_height   = GRID_ROW_HEIGHT,
                          units        = "cm",
                          x_dot        = GRID_X_DOT,
                          gap_dot_num  = GRID_GAP_DN,
                          gap_num_name = GRID_GAP_NN,
                          dot_alpha    = GRID_DOT_ALPHA,
                          isf_note     = ISF_NOTE_MAIN,
                          entry_gap    = GRID_ENTRY_GAP,
                          title_gap    = GRID_TITLE_GAP,
                          section_gap  = GRID_SECTION_GAP,
                          y_npc        = GRID_LEG_Y_NPC)
  
  leg_gg <- ggdraw() + draw_grob(leg) +
    theme(plot.background = element_rect(fill = "transparent", color = NA))
  
  leg_w <- legend_width_cm(sub_df_leg,
                           text_size    = GRID_TEXT_SIZE,
                           title_size   = GRID_LEG_TITLE_SIZE,
                           x_dot        = GRID_X_DOT,
                           gap_dot_num  = GRID_GAP_DN,
                           gap_num_name = GRID_GAP_NN,
                           pad_cm       = GRID_LEG_PAD_CM)
  
  cell <- if (is.null(break_range)) {
    plot_grid(p, leg_gg, nrow = 1,
              rel_widths = c(GRID_PANEL_W_CM, leg_w), align = "h", axis = "b")
  } else {
    plot_grid(p, leg_gg, nrow = 1,
              rel_widths = c(GRID_PANEL_W_CM, leg_w))
  }
  
  list(plot = cell, width = GRID_PANEL_W_CM + leg_w)
}

# ==============================================================================
# 13. Build and save the main grid -----
# ==============================================================================
break_groups <- names(xbreak_settings)
solo_groups  <- "PC(16:0_9:0)"
main_groups  <- setdiff(unique(For_mass_defect$Group3), c(break_groups, solo_groups))

panel_list <- lapply(main_groups, build_panel)
cells      <- lapply(panel_list, `[[`, "plot")
cell_w     <- vapply(panel_list, `[[`, numeric(1), "width")

n_cols  <- 4
row_idx <- split(seq_along(cells), ceiling(seq_along(cells) / n_cols))

row_plots <- lapply(row_idx, function(i) {
  plot_grid(plotlist = cells[i], nrow = 1, rel_widths = cell_w[i])
})

final_plot <- if (length(row_plots) == 1) {
  row_plots[[1]]
} else {
  plot_grid(plotlist = row_plots, ncol = 1)
}

final_plot

# ==============================================================================
# 13b. Legend bottom alignment -----
# ==============================================================================
panel_bottom_pad <- function(p) {
  gt      <- ggplotGrob(p)
  panel_b <- max(gt$layout$b[grepl("^panel", gt$layout$name)])
  sum(convertHeight(gt$heights[(panel_b + 1):length(gt$heights)], "cm", valueOnly = TRUE))
}

# ==============================================================================
# 14. PC(16:0_18:1) with axis break -----
# ==============================================================================
grp         <- "PC(16:0_18:1)"
break_range <- xbreak_settings[[grp]]$gap
scale_ratio <- xbreak_settings[[grp]]$scale_ratio
xlim_range  <- xlim_settings[[grp]]

sub_df <- For_mass_defect %>%
  filter(Group3 == grp) %>%
  as.data.frame()

sub_df_before_break <- sub_df %>% filter(Average.Rt.min. <= break_range[1])
sub_df_after_break  <- sub_df %>% filter(Average.Rt.min. >= break_range[2])

x_compensation <- if (is.numeric(scale_ratio)) scale_ratio else 1

sub_df_before_break <- position_labels_by_kmd(sub_df_before_break,
                                              relative_label_width  = 0.08,
                                              relative_spacing      = 0.2,
                                              relative_y_offset     = 0.05,
                                              cluster_gap_threshold = 0.15)

sub_df_after_break <- position_labels_by_kmd(sub_df_after_break,
                                             relative_label_width  = 0.4,
                                             relative_spacing      = 0.12 * x_compensation,
                                             relative_y_offset     = 0.8,
                                             cluster_gap_threshold = 0.15)

sub_df <- bind_rows(sub_df_before_break, sub_df_after_break)

sub_df <- sub_df %>%
  mutate(precursor_label = ifelse(color_group == "Non-oxidized",
                                  str_extract(Metabolite.curated, "\\([a-z]\\)$"),
                                  NA_character_))

carbondb_values  <- unique(sub_df$CarbonDB_for_shape[sub_df$color_group != "ISF"])
carbondb_values  <- carbondb_values[!is.na(carbondb_values)]
carbondb_sorted  <- c(intersect("matched", carbondb_values),
                      setdiff(carbondb_values, "matched") %>% sort())
shape_vals_panel <- setNames(shape_values[seq_along(carbondb_sorted)], carbondb_sorted)

sub_df <- sub_df %>%
  mutate(CarbonDB_for_shape = ifelse(color_group == "ISF", "ISF_marker", CarbonDB_for_shape))
shape_vals_panel["ISF_marker"] <- ISF_SHAPE

colors_used          <- unique(sub_df$color_group)
fill_colors_panel    <- my_fill_colors[colors_used]
outline_colors_panel <- my_outline_colors[colors_used]

n_precursors   <- sum(sub_df$oxidation_status == "Non-oxidized")
precursor_word <- ifelse(n_precursors > 1, "precursors", "precursor")
panel_title    <- paste0("Ox", grp, " and ", precursor_word)

p_nobreak <- ggplot(sub_df,
                    aes(x = Average.Rt.min., y = KMD,
                        fill = color_group, color = color_group, shape = CarbonDB_for_shape)) +
  geom_segment(
    data = subset(sub_df, !(color_group == "Non-oxidized" & is.na(precursor_label)) & color_group != "ISF"),
    aes(x = Average.Rt.min., y = KMD, xend = label_x, yend = label_y),
    color = "grey40", size = 0.3) +
  geom_point(size = SOLO_PT_SIZE, alpha = SOLO_DOT_ALPHA, stroke = 0.5) +
  geom_text(data = subset(sub_df, color_group != "Non-oxidized" & !is.na(point_label)),
            aes(x = label_x, y = label_y, label = point_label),
            color = "black", size = SOLO_PT_LABEL_SIZE / .pt, hjust = 0.5, vjust = -0.2) +
  geom_text(data = subset(sub_df, color_group == "Non-oxidized" & !is.na(precursor_label)),
            aes(x = label_x, y = label_y, label = precursor_label),
            color = "black", size = SOLO_PT_LABEL_SIZE / .pt, hjust = 0.5, vjust = -0.2) +
  scale_fill_manual(values = fill_colors_panel) +
  scale_color_manual(values = outline_colors_panel) +
  scale_shape_manual(values = shape_vals_panel) +
  scale_x_continuous(labels = function(x) sprintf("%.2f", x), limits = xlim_range) +
  geom_segment(data = subset(sub_df, color_group == "ISF"),
               aes(x    = Average.Rt.min. - ISF_LINE_WIDTH_MIN / 2,
                   xend = Average.Rt.min. + ISF_LINE_WIDTH_MIN / 2,
                   y    = KMD + ISF_LINE_Y_OFFSET,
                   yend = KMD + ISF_LINE_Y_OFFSET),
               inherit.aes = FALSE,
               color = "red", linetype = ISF_LINE_TYPE, linewidth = ISF_LINE_WIDTH) +
  scale_y_continuous(expand = expansion(mult = 0.25),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  labs(title = panel_title, x = "Retention time (min)", y = "O-KMD") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title           = element_text(face = "bold", size = SOLO_TITLE_SIZE, color = "black",
                                        hjust = 0.65, vjust = -3.5),
    axis.title.x         = element_text(face = "bold", size = SOLO_AXIS_TITLE_SIZE, color = "black",
                                        hjust = 0.6, margin = margin(t = -1, unit = "pt")),
    axis.title.y         = element_text(face = "bold", size = SOLO_AXIS_TITLE_SIZE, color = "black",
                                        hjust = 0.5, margin = margin(r = SOLO_YTITLE_MARGIN, unit = "pt")),
    axis.text.x          = element_text(size = SOLO_AXIS_TEXT_SIZE, color = "black",
                                        margin = margin(t = 8, unit = "pt")),
    axis.text.x.top      = element_blank(),
    axis.ticks.x.top     = element_blank(),
    axis.text.y          = element_text(size = SOLO_AXIS_TEXT_SIZE, color = "black",
                                        margin = margin(r = 8, unit = "pt")),
    legend.position      = "none",
    panel.grid.major     = element_blank(),
    panel.grid.minor     = element_blank(),
    panel.border         = element_blank(),
    axis.line            = element_line(color = "black", linewidth = 0.2),
    axis.line.y.right    = element_line(color = "black", linewidth = 0.2),
    axis.ticks           = element_line(color = "black", linewidth = 0.3),
    axis.ticks.y.right   = element_blank(),
    axis.ticks.length.x  = unit(-SOLO_TICK_LENGTH, "cm"),
    axis.ticks.length.y  = unit(-SOLO_TICK_LENGTH, "cm"),
    plot.margin          = margin(t = 4, r = 4, b = 4, l = 2),
    plot.background      = element_rect(fill = "transparent", color = NA),
    panel.background     = element_rect(fill = "transparent", color = NA)
  )

leg_pad_cm <- panel_bottom_pad(p_nobreak)

p <- p_nobreak +
  scale_x_break(break_range, scales = scale_ratio,
                ticklabels = c(14.15, 20.15))

sub_df_leg <- sub_df %>%
  mutate(CarbonDB_for_shape = ifelse(color_group == "ISF", NA_character_, CarbonDB_for_shape)) %>%
  mutate(
    Metabolite.curated = ifelse(color_group == "Non-oxidized",
                                str_remove(Metabolite.curated, "\\([a-z]\\)$"),
                                Metabolite.curated),
    point_label = ifelse(color_group == "Non-oxidized", NA_integer_, point_label)
  ) %>%
  group_by(Metabolite.curated) %>%
  mutate(
    Metabolite_leg = ifelse(
      color_group == "Non-oxidized" & any(!is.na(precursor_label)),
      paste0(Metabolite.curated,
             paste(sort(unique(precursor_label[!is.na(precursor_label)])), collapse = "/")),
      Metabolite.curated
    )
  ) %>%
  ungroup() %>%
  mutate(Metabolite.curated = Metabolite_leg) %>%
  select(-Metabolite_leg)

leg <- make_legend_grob(sub_df_leg,
                        dot_size     = SOLO_DOT_SIZE,
                        text_size    = SOLO_TEXT_SIZE,
                        title_size   = SOLO_LEG_TITLE_SIZE,
                        row_height   = SOLO_ROW_HEIGHT,
                        gap_dot_num  = SOLO_GAP_DN,
                        gap_num_name = SOLO_GAP_NN,
                        dot_alpha    = SOLO_DOT_ALPHA,
                        isf_note     = ISF_NOTE_BREAK,
                        title_nonox  = "Precursor lipids",
                        entry_gap    = SOLO_ENTRY_GAP,
                        title_gap    = SOLO_TITLE_GAP,
                        section_gap  = SOLO_SECTION_GAP,
                        y_npc        = SOLO_LEG_Y_NPC)

LEG_EXTRA_PAD_CM <- 0.18

leg_gg <- ggdraw() + draw_grob(leg) +
  theme(plot.margin     = margin(t = 0, r = 0, b = leg_pad_cm + LEG_EXTRA_PAD_CM, l = -0.38, unit = "cm"),
        plot.background = element_rect(fill = "transparent", color = NA))

p_grob <- ggdraw() +
  draw_grob(make_grob_transparent(grid.grabExpr(print(p)))) +
  theme(plot.background = element_rect(fill = "transparent", color = NA))

final_PC_18_1 <- plot_grid(p_grob, leg_gg, nrow = 1, rel_widths = c(1, SOLO_LEG_WIDTH))
final_PC_18_1

# ==============================================================================
# 15. PC(16:0_16:1) with axis break -----
# ==============================================================================
grp         <- "PC(16:0_16:1)"
break_range <- xbreak_settings[[grp]]$gap
scale_ratio <- xbreak_settings[[grp]]$scale_ratio
xlim_range  <- xlim_settings[[grp]]

sub_df <- For_mass_defect %>%
  filter(Group3 == grp) %>%
  as.data.frame()

sub_df_before_break <- sub_df %>% filter(Average.Rt.min. <= break_range[1])
sub_df_after_break  <- sub_df %>% filter(Average.Rt.min. >= break_range[2])

x_compensation <- if (is.numeric(scale_ratio)) scale_ratio else 1

sub_df_before_break <- position_labels_by_kmd(sub_df_before_break,
                                              relative_label_width  = 0.05,
                                              relative_spacing      = 0.08,
                                              relative_y_offset     = 0.05,
                                              cluster_gap_threshold = 0.15)

sub_df_after_break <- position_labels_by_kmd(sub_df_after_break,
                                             relative_label_width  = 0.05 * x_compensation,
                                             relative_spacing      = 0.07 * x_compensation,
                                             relative_y_offset     = 0.15,
                                             cluster_gap_threshold = 0.15)

sub_df <- bind_rows(sub_df_before_break, sub_df_after_break)

sub_df <- sub_df %>%
  mutate(precursor_label = ifelse(color_group == "Non-oxidized",
                                  str_extract(Metabolite.curated, "\\([a-z]\\)$"),
                                  NA_character_))

carbondb_values  <- unique(sub_df$CarbonDB_for_shape[sub_df$color_group != "ISF"])
carbondb_values  <- carbondb_values[!is.na(carbondb_values)]
carbondb_sorted  <- c(intersect("matched", carbondb_values),
                      setdiff(carbondb_values, "matched") %>% sort())
shape_vals_panel <- setNames(shape_values[seq_along(carbondb_sorted)], carbondb_sorted)

sub_df <- sub_df %>%
  mutate(CarbonDB_for_shape = ifelse(color_group == "ISF", "ISF_marker", CarbonDB_for_shape))
shape_vals_panel["ISF_marker"] <- ISF_SHAPE

colors_used          <- unique(sub_df$color_group)
fill_colors_panel    <- my_fill_colors[colors_used]
outline_colors_panel <- my_outline_colors[colors_used]

n_precursors   <- sum(sub_df$oxidation_status == "Non-oxidized")
precursor_word <- ifelse(n_precursors > 1, "precursors", "precursor")
panel_title    <- paste0("Ox", grp, " and ", precursor_word)

p_nobreak <- ggplot(sub_df,
                    aes(x = Average.Rt.min., y = KMD,
                        fill = color_group, color = color_group, shape = CarbonDB_for_shape)) +
  geom_segment(
    data = subset(sub_df, !(color_group == "Non-oxidized" & is.na(precursor_label)) & color_group != "ISF"),
    aes(x = Average.Rt.min., y = KMD, xend = label_x, yend = label_y),
    color = "grey40", size = 0.3) +
  geom_point(size = SOLO_PT_SIZE, alpha = SOLO_DOT_ALPHA, stroke = 0.5) +
  geom_text(data = subset(sub_df, color_group != "Non-oxidized" & !is.na(point_label)),
            aes(x = label_x, y = label_y, label = point_label),
            color = "black", size = SOLO_PT_LABEL_SIZE / .pt, hjust = 0.5, vjust = -0.2) +
  geom_text(data = subset(sub_df, color_group == "Non-oxidized" & !is.na(precursor_label)),
            aes(x = label_x, y = label_y, label = precursor_label),
            color = "black", size = SOLO_PT_LABEL_SIZE / .pt, hjust = 0.5, vjust = -0.2) +
  scale_fill_manual(values = fill_colors_panel) +
  scale_color_manual(values = outline_colors_panel) +
  scale_shape_manual(values = shape_vals_panel) +
  scale_x_continuous(labels = function(x) sprintf("%.2f", x), limits = xlim_range) +
  geom_segment(data = subset(sub_df, color_group == "ISF"),
               aes(x    = Average.Rt.min. - ISF_LINE_WIDTH_MIN / 2,
                   xend = Average.Rt.min. + ISF_LINE_WIDTH_MIN / 2,
                   y    = KMD + ISF_LINE_Y_OFFSET,
                   yend = KMD + ISF_LINE_Y_OFFSET),
               inherit.aes = FALSE,
               color = "red", linetype = ISF_LINE_TYPE, linewidth = ISF_LINE_WIDTH) +
  scale_y_continuous(expand = expansion(mult = 0.25),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  labs(title = panel_title, x = "Retention time (min)", y = "O-KMD") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title           = element_text(face = "bold", size = SOLO_TITLE_SIZE, color = "black",
                                        hjust = 0.65, vjust = -3.5),
    axis.title.x         = element_text(face = "bold", size = SOLO_AXIS_TITLE_SIZE, color = "black",
                                        hjust = 0.6, margin = margin(t = -1, unit = "pt")),
    axis.title.y         = element_text(face = "bold", size = SOLO_AXIS_TITLE_SIZE, color = "black",
                                        hjust = 0.5, margin = margin(r = SOLO_YTITLE_MARGIN, unit = "pt")),
    axis.text.x          = element_text(size = SOLO_AXIS_TEXT_SIZE, color = "black",
                                        margin = margin(t = 8, unit = "pt")),
    axis.text.x.top      = element_blank(),
    axis.ticks.x.top     = element_blank(),
    axis.text.y          = element_text(size = SOLO_AXIS_TEXT_SIZE, color = "black",
                                        margin = margin(r = 8, unit = "pt")),
    legend.position      = "none",
    panel.grid.major     = element_blank(),
    panel.grid.minor     = element_blank(),
    panel.border         = element_blank(),
    axis.line            = element_line(color = "black", linewidth = 0.2),
    axis.line.y.right    = element_line(color = "black", linewidth = 0.2),
    axis.ticks           = element_line(color = "black", linewidth = 0.3),
    axis.ticks.y.right   = element_blank(),
    axis.ticks.length.x  = unit(-SOLO_TICK_LENGTH, "cm"),
    axis.ticks.length.y  = unit(-SOLO_TICK_LENGTH, "cm"),
    plot.margin          = margin(t = 4, r = 4, b = 4, l = 2),
    plot.background      = element_rect(fill = "transparent", color = NA),
    panel.background     = element_rect(fill = "transparent", color = NA)
  )

leg_pad_cm <- panel_bottom_pad(p_nobreak)

p <- p_nobreak +
  scale_x_break(break_range, scales = scale_ratio,
                ticklabels = c(11.45, 13.45))

sub_df_leg <- sub_df %>%
  mutate(CarbonDB_for_shape = ifelse(color_group == "ISF", NA_character_, CarbonDB_for_shape)) %>%
  mutate(
    Metabolite.curated = ifelse(color_group == "Non-oxidized",
                                str_remove(Metabolite.curated, "\\([a-z]\\)$"),
                                Metabolite.curated),
    point_label = ifelse(color_group == "Non-oxidized", NA_integer_, point_label)
  ) %>%
  group_by(Metabolite.curated) %>%
  mutate(
    Metabolite_leg = ifelse(
      color_group == "Non-oxidized" & any(!is.na(precursor_label)),
      paste0(Metabolite.curated,
             paste(sort(unique(precursor_label[!is.na(precursor_label)])), collapse = "/")),
      Metabolite.curated
    )
  ) %>%
  ungroup() %>%
  mutate(Metabolite.curated = Metabolite_leg) %>%
  select(-Metabolite_leg)

leg <- make_legend_grob(sub_df_leg,
                        dot_size     = SOLO_DOT_SIZE,
                        text_size    = SOLO_TEXT_SIZE,
                        title_size   = SOLO_LEG_TITLE_SIZE,
                        row_height   = SOLO_ROW_HEIGHT,
                        gap_dot_num  = SOLO_GAP_DN,
                        gap_num_name = SOLO_GAP_NN,
                        dot_alpha    = SOLO_DOT_ALPHA,
                        isf_note     = ISF_NOTE_BREAK,
                        title_nonox  = "Precursor lipids",
                        entry_gap    = SOLO_ENTRY_GAP,
                        title_gap    = SOLO_TITLE_GAP,
                        section_gap  = SOLO_SECTION_GAP,
                        y_npc        = SOLO_LEG_Y_NPC)

LEG_EXTRA_PAD_CM <- 0.18

leg_gg <- ggdraw() + draw_grob(leg) +
  theme(plot.margin     = margin(t = 0, r = 0, b = leg_pad_cm + LEG_EXTRA_PAD_CM, l = -0.38, unit = "cm"),
        plot.background = element_rect(fill = "transparent", color = NA))

p_grob <- ggdraw() +
  draw_grob(make_grob_transparent(grid.grabExpr(print(p)))) +
  theme(plot.background = element_rect(fill = "transparent", color = NA))

final_PC_16_1 <- plot_grid(p_grob, leg_gg, nrow = 1, rel_widths = c(1, SOLO_LEG_WIDTH))
final_PC_16_1

# ==============================================================================
# 16. PC(16:0_9:0) standalone panel, no axis break -----
# ==============================================================================
grp        <- "PC(16:0_9:0)"
xlim_range <- xlim_settings[[grp]]

SOLO_LEG_LEFT_PAD_CM  <- -0.1
SOLO_LEG_EXTRA_PAD_CM <- -0
PC90_YTITLE_MARGIN    <- 0

sub_df <- For_mass_defect %>%
  filter(Group3 == grp) %>%
  as.data.frame()

sub_df <- position_labels_by_kmd(sub_df,
                                 relative_label_width  = 0.05,
                                 relative_spacing      = 0.08,
                                 relative_y_offset     = 0.04,
                                 cluster_gap_threshold = 0.15)

sub_df <- sub_df %>%
  mutate(precursor_label = ifelse(color_group == "Non-oxidized",
                                  str_extract(Metabolite.curated, "\\([a-z]\\)$"),
                                  NA_character_))

carbondb_values  <- unique(sub_df$CarbonDB_for_shape[sub_df$color_group != "ISF"])
carbondb_values  <- carbondb_values[!is.na(carbondb_values)]
carbondb_sorted  <- c(intersect("matched", carbondb_values),
                      setdiff(carbondb_values, "matched") %>% sort())
shape_vals_panel <- setNames(shape_values[seq_along(carbondb_sorted)], carbondb_sorted)

sub_df <- sub_df %>%
  mutate(CarbonDB_for_shape = ifelse(color_group == "ISF", "ISF_marker", CarbonDB_for_shape))
shape_vals_panel["ISF_marker"] <- ISF_SHAPE

colors_used          <- unique(sub_df$color_group)
fill_colors_panel    <- my_fill_colors[colors_used]
outline_colors_panel <- my_outline_colors[colors_used]

n_precursors   <- sum(sub_df$oxidation_status == "Non-oxidized")
precursor_word <- ifelse(n_precursors > 1, "precursors", "precursor")
panel_title    <- paste0("Ox", grp, " and ", precursor_word)

p_solo <- ggplot(sub_df,
                 aes(x = Average.Rt.min., y = KMD,
                     fill = color_group, color = color_group, shape = CarbonDB_for_shape)) +
  geom_segment(
    data = subset(sub_df, !(color_group == "Non-oxidized" & is.na(precursor_label)) & color_group != "ISF"),
    aes(x = Average.Rt.min., y = KMD, xend = label_x, yend = label_y),
    color = "grey40", size = 0.3) +
  geom_point(size = SOLO_PT_SIZE, alpha = SOLO_DOT_ALPHA, stroke = 0.5) +
  geom_text(data = subset(sub_df, color_group != "Non-oxidized" & !is.na(point_label)),
            aes(x = label_x, y = label_y, label = point_label),
            color = "black", size = SOLO_PT_LABEL_SIZE / .pt, hjust = 0.5, vjust = -0.2) +
  geom_text(data = subset(sub_df, color_group == "Non-oxidized" & !is.na(precursor_label)),
            aes(x = label_x, y = label_y, label = precursor_label),
            color = "black", size = SOLO_PT_LABEL_SIZE / .pt, hjust = 0.5, vjust = -0.2) +
  scale_fill_manual(values = fill_colors_panel) +
  scale_color_manual(values = outline_colors_panel) +
  scale_shape_manual(values = shape_vals_panel) +
  scale_x_continuous(labels = function(x) sprintf("%.2f", x), limits = xlim_range) +
  scale_y_continuous(expand = expansion(mult = 0.25)) +
  labs(title = panel_title, x = "Retention time (min)", y = "O-KMD") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title          = element_text(face = "bold", size = SOLO_TITLE_SIZE, color = "black", hjust = 0.5,
                                       margin = margin(b = 2, unit = "pt")),
    axis.title.x        = element_text(face = "bold", size = SOLO_AXIS_TITLE_SIZE, color = "black"),
    axis.title.y        = element_text(face = "bold", size = SOLO_AXIS_TITLE_SIZE, color = "black",
                                       margin = margin(r = PC90_YTITLE_MARGIN, unit = "pt")),
    axis.text.x         = element_text(size = SOLO_AXIS_TEXT_SIZE, color = "black",
                                       margin = margin(t = 8, unit = "pt")),
    axis.text.y         = element_text(size = SOLO_AXIS_TEXT_SIZE, color = "black",
                                       margin = margin(r = 8, unit = "pt")),
    legend.position     = "none",
    panel.grid.major    = element_blank(),
    panel.grid.minor    = element_blank(),
    panel.border        = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks          = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length.x = unit(-SOLO_TICK_LENGTH, "cm"),
    axis.ticks.length.y = unit(-SOLO_TICK_LENGTH, "cm"),
    plot.margin         = margin(t = 4, r = 4, b = 4, l = 2),
    plot.background     = element_rect(fill = "transparent", color = NA),
    panel.background    = element_rect(fill = "transparent", color = NA)
  )

leg_pad_cm <- panel_bottom_pad(p_solo)

sub_df_leg <- sub_df %>%
  mutate(CarbonDB_for_shape = ifelse(color_group == "ISF", NA_character_, CarbonDB_for_shape)) %>%
  mutate(
    Metabolite.curated = ifelse(color_group == "Non-oxidized",
                                str_remove(Metabolite.curated, "\\([a-z]\\)$"),
                                Metabolite.curated),
    point_label = ifelse(color_group == "Non-oxidized", NA_integer_, point_label)
  ) %>%
  group_by(Metabolite.curated) %>%
  mutate(
    Metabolite_leg = ifelse(
      color_group == "Non-oxidized" & any(!is.na(precursor_label)),
      paste0(Metabolite.curated,
             paste(sort(unique(precursor_label[!is.na(precursor_label)])), collapse = "/")),
      Metabolite.curated
    )
  ) %>%
  ungroup() %>%
  mutate(Metabolite.curated = Metabolite_leg) %>%
  select(-Metabolite_leg)

leg <- make_legend_grob(sub_df_leg,
                        dot_size     = SOLO_DOT_SIZE,
                        text_size    = SOLO_TEXT_SIZE,
                        title_size   = SOLO_LEG_TITLE_SIZE,
                        row_height   = SOLO_ROW_HEIGHT,
                        gap_dot_num  = SOLO_GAP_DN,
                        gap_num_name = SOLO_GAP_NN,
                        dot_alpha    = SOLO_DOT_ALPHA,
                        isf_note     = ISF_NOTE_MAIN,
                        title_nonox  = "Precursor lipids",
                        entry_gap    = SOLO_ENTRY_GAP,
                        title_gap    = SOLO_TITLE_GAP,
                        section_gap  = SOLO_SECTION_GAP,
                        y_npc        = SOLO_LEG_Y_NPC)

leg_gg <- ggdraw() + draw_grob(leg) +
  theme(plot.margin     = margin(t = 0, r = 0,
                                 b = leg_pad_cm + SOLO_LEG_EXTRA_PAD_CM,
                                 l = SOLO_LEG_LEFT_PAD_CM, unit = "cm"),
        plot.background = element_rect(fill = "transparent", color = NA))

final_PC_9_0 <- plot_grid(p_solo, leg_gg, nrow = 1, rel_widths = c(1, SOLO_LEG_WIDTH))
final_PC_9_0

# ==============================================================================
# 17. Save -----
# ==============================================================================
ggsave("./Figure/O-Kendrick_mass_defect_plots.pdf",
       plot   = final_plot,
       width  = 1230*0.95/96,
       height = 240*0.9/96,
       units  = "in",
       bg     = "transparent")

ggsave("./Figure/O-KMD_PC_16_0_18_1_break.pdf",
       plot   = final_PC_18_1,
       width  = 500*1.25/96,
       height = 280*1/96,
       units  = "in",
       bg     = "transparent")

ggsave("./Figure/O-KMD_PC_16_0_16_1_break.pdf",
       plot   = final_PC_16_1,
       width  = 480*1.2/96,
       height = 280*1/96,
       units  = "in",
       bg     = "transparent")

ggsave("./Figure/O-KMD_PC_16_0_9_0.pdf",
       plot   = final_PC_9_0,
       width  = 480*1.2*0.95/96,
       height = 275*1*0.95/96,
       units  = "in",
       bg     = "transparent")