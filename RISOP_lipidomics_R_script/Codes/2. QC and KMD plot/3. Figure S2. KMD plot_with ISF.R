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

# 1.1 Non-oxidized lipids
NonOx_df <- read.csv('./Input/NonOx_df.csv', stringsAsFactors = F, skip = 0, header = T)

# 1.2 Oxidized lipids
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

ISF_IDs <- c("38015", "39402", "30767", "39995", "41440", "39887", "42220")
Ox_df_curated <- Ox_df_curated %>%
  mutate(is_ISF = Alignment.ID %in% ISF_IDs)

# 1.3 Library file
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
base_nominal <- 1.000000
base_exact   <- 1.007825032
scale_factor <- base_nominal / base_exact

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
                             row_height   = 0.4,
                             gap_dot_num  = 0.08,
                             gap_num_name = 0.06,
                             dot_alpha    = 0.8) {
  
  carbondb_values <- unique(panel_data$CarbonDB_for_shape[panel_data$color_group != "ISF"])
  carbondb_values <- carbondb_values[!is.na(carbondb_values)]
  
  carbondb_sorted <- c(
    intersect("matched", carbondb_values),
    setdiff(carbondb_values, "matched") %>% sort()
  )
  
  pch_lookup <- setNames(shape_values[seq_along(carbondb_sorted)], carbondb_sorted)
  
  # ISF lipids sorted to the bottom
  legend_df <- panel_data %>%
    arrange(color_group == "ISF", point_label) %>%
    select(point_label, color_group, Metabolite.curated, CarbonDB_for_shape) %>%
    distinct()
  
  # ISF rows occupy 2 units of row_height; non-ISF rows occupy 1 unit
  row_weights      <- ifelse(legend_df$color_group == "ISF", 2, 1)
  total_rows       <- sum(row_weights)
  total_height_cm  <- total_rows * row_height
  
  n <- nrow(legend_df)
  
  x_dot  <- 0.05
  x_num  <- x_dot + gap_dot_num
  x_name <- x_num + gap_num_name
  
  grob_list <- list()
  
  # Track cumulative y position using row_weights
  cumulative_units <- 0
  
  for (i in seq_len(n)) {
    row_i      <- legend_df[i, ]
    is_isf_row <- row_i$color_group == "ISF"
    w          <- row_weights[i]
    
    row_center_units <- cumulative_units + w / 2
    y_pos            <- 1 - row_center_units / total_rows
    
    cumulative_units <- cumulative_units + w
    
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
    
    if (is_isf_row) {
      # Dot and lipid name in upper half; sub-label in lower half of the 2-unit space
      upper_y <- 1 - (cumulative_units - w + 0.5) / total_rows
      lower_y <- 1 - (cumulative_units - 0.5)     / total_rows
      
      grob_list[[length(grob_list) + 1]] <- pointsGrob(
        x    = unit(x_dot, "npc"),
        y    = unit(upper_y, "npc"),
        pch  = pch_val,
        size = unit(dot_size, "pt"),
        gp   = gpar(fill = fill_col, col = outline_col, lwd = 0.8)
      )
      grob_list[[length(grob_list) + 1]] <- textGrob(
        label = as.character(row_i$point_label),
        x     = unit(x_num, "npc"),
        y     = unit(upper_y, "npc"),
        just  = "left",
        gp    = gpar(fontsize = text_size, col = "black")
      )
      grob_list[[length(grob_list) + 1]] <- textGrob(
        label = row_i$Metabolite.curated,
        x     = unit(x_name, "npc"),
        y     = unit(upper_y, "npc"),
        just  = "left",
        gp    = gpar(fontsize = text_size, col = "black")
      )
      grob_list[[length(grob_list) + 1]] <- textGrob(
        label = "(potential ISF product, removed)",
        x     = unit(x_name, "npc"),
        y     = unit(lower_y, "npc"),
        just  = "left",
        gp    = gpar(fontsize = text_size, col = "black")
      )
      
    } else {
      grob_list[[length(grob_list) + 1]] <- pointsGrob(
        x    = unit(x_dot, "npc"),
        y    = unit(y_pos, "npc"),
        pch  = pch_val,
        size = unit(dot_size, "pt"),
        gp   = gpar(fill = fill_col, col = outline_col, lwd = 0.8)
      )
      grob_list[[length(grob_list) + 1]] <- textGrob(
        label = as.character(row_i$point_label),
        x     = unit(x_num, "npc"),
        y     = unit(y_pos, "npc"),
        just  = "left",
        gp    = gpar(fontsize = text_size, col = "black")
      )
      grob_list[[length(grob_list) + 1]] <- textGrob(
        label = row_i$Metabolite.curated,
        x     = unit(x_name, "npc"),
        y     = unit(y_pos, "npc"),
        just  = "left",
        gp    = gpar(fontsize = text_size, col = "black")
      )
    }
  }
  
  gTree(
    children = do.call(gList, grob_list),
    vp = viewport(
      width  = unit(1, "npc"),
      height = unit(total_height_cm, "cm"),
      y      = unit(0, "npc"),
      just   = c("center", "bottom")
    )
  )
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
  mutate(KMD_rounded = round(KMD / 0.05) * 0.05) %>%
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
  mutate(
    n_unique_kmd = length(unique(KMD_row_id[!is.na(KMD_row_id)]))
  ) %>%
  ungroup() %>%
  group_by(Group3, KMD_row_id) %>%
  mutate(
    min_RT_in_row = ifelse(oxidation_status == "Oxidized",
                           min(Average.Rt.min.),
                           NA)
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
  mutate(
    min_RT_per_kmd_group = min(Average.Rt.min.)
  ) %>%
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
  mutate(
    Group3_CarbonDB = str_extract(Group3, "(?<=\\()[0-9:_]+(?=\\))")
  ) %>%
  group_by(Group3) %>%
  mutate(
    CarbonDB_for_shape = case_when(
      color_group == "ISF"            ~ NA_character_,
      CarbonDB == Group3_CarbonDB     ~ "matched",
      TRUE                            ~ CarbonDB
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
            seq(center_x - half_span,
                center_x + half_span,
                length.out = n)
          }
        }
      }
    ) %>%
    ungroup() %>%
    select(point_label, label_x, label_y)
  
  sub_df <- sub_df %>%
    left_join(label_positions, by = "point_label")
  
  return(sub_df)
}

# ==============================================================================
# 11. Global plot tuning parameters -----
# ==============================================================================
ALL_DOT_SIZE    <- 6
ALL_TEXT_SIZE   <- 6
ALL_ROW_HEIGHT  <- 0.32
ALL_GAP_DN      <- 0.05
ALL_GAP_NN      <- 0.10
ALL_DOT_ALPHA   <- 0.8
ALL_LEG_WIDTH   <- 0.55
ALL_TICK_LENGTH <- 0.05

# ==============================================================================
# 12. Build one panel per lipid group -----
# ==============================================================================
all_groups <- unique(For_mass_defect$Group3)
panel_list <- lapply(all_groups, function(grp) {
  
  sub_df <- For_mass_defect %>%
    filter(Group3 == grp) %>%
    as.data.frame()
  
  sub_df <- position_labels_by_kmd(sub_df,
                                   relative_label_width  = 0.05,
                                   relative_spacing      = 0.08,
                                   relative_y_offset     = 0.06,
                                   cluster_gap_threshold = 0.15)
  
  carbondb_values <- unique(sub_df$CarbonDB_for_shape[sub_df$color_group != "ISF"])
  carbondb_values <- carbondb_values[!is.na(carbondb_values)]
  
  carbondb_sorted <- c(
    intersect("matched", carbondb_values),
    setdiff(carbondb_values, "matched") %>% sort()
  )
  
  shape_vals_panel <- setNames(shape_values[seq_along(carbondb_sorted)], carbondb_sorted)
  
  sub_df <- sub_df %>%
    mutate(CarbonDB_for_shape = ifelse(color_group == "ISF", "ISF_marker", CarbonDB_for_shape))
  
  shape_vals_panel["ISF_marker"] <- ISF_SHAPE
  
  colors_used          <- unique(sub_df$color_group)
  fill_colors_panel    <- my_fill_colors[colors_used]
  outline_colors_panel <- my_outline_colors[colors_used]
  
  n_precursors   <- sum(sub_df$oxidation_status == "Non-oxidized")
  precursor_word <- ifelse(n_precursors > 1, "precursors", "precursor")
  
  p <- ggplot(sub_df,
              aes(x = Average.Rt.min., y = KMD,
                  fill = color_group, color = color_group, shape = CarbonDB_for_shape)) +
    geom_segment(aes(x = Average.Rt.min., y = KMD,
                     xend = label_x, yend = label_y),
                 color = "grey40", size = 0.3) +
    geom_point(size = 2, alpha = ALL_DOT_ALPHA, stroke = 0.5) +
    geom_text(aes(x = label_x, y = label_y, label = point_label),
              color = "black", size = 2, hjust = 0.5, vjust = -0.2) +
    scale_fill_manual(values = fill_colors_panel) +
    scale_color_manual(values = outline_colors_panel) +
    scale_shape_manual(values = shape_vals_panel) +
    scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
    scale_y_continuous(expand = expansion(mult = 0.25)) +
    labs(title = paste0("Ox", grp, " and ", precursor_word),
         x = "Retention time (min)", y = "H-KMD") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title          = element_text(face = "bold", size = 8, color = "black", hjust = 0.5,
                                         margin = margin(b = -0.0, unit = "pt")),
      axis.title.x        = element_text(face = "bold", size = 8, color = "black"),
      axis.title.y        = element_text(face = "bold", size = 8, color = "black"),
      axis.text.x         = element_text(size = 7, color = "black",
                                         margin = margin(t = 8, unit = "pt")),
      axis.text.y         = element_text(size = 7, color = "black",
                                         margin = margin(r = 8, unit = "pt")),
      legend.position     = "none",
      panel.grid.major    = element_blank(),
      panel.grid.minor    = element_blank(),
      panel.border        = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.ticks          = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length.x = unit(-ALL_TICK_LENGTH, "cm"),
      axis.ticks.length.y = unit(-ALL_TICK_LENGTH, "cm"),
      plot.margin         = margin(t = 4, r = 4, b = 4, l = 2)
    )
  
  sub_df_leg <- sub_df %>%
    mutate(CarbonDB_for_shape = ifelse(color_group == "ISF", NA_character_, CarbonDB_for_shape))
  
  leg    <- make_legend_grob(sub_df_leg,
                             dot_size     = ALL_DOT_SIZE,
                             text_size    = ALL_TEXT_SIZE,
                             row_height   = ALL_ROW_HEIGHT,
                             gap_dot_num  = ALL_GAP_DN,
                             gap_num_name = ALL_GAP_NN,
                             dot_alpha    = ALL_DOT_ALPHA)
  leg_gg <- ggdraw() + draw_grob(leg)
  
  plot_grid(p, leg_gg, nrow = 1, rel_widths = c(1, ALL_LEG_WIDTH), align = "h", axis = "b")
})

# ==============================================================================
# 13. Assemble and save final figure -----
# ==============================================================================
n_cols     <- 3
final_plot <- do.call(plot_grid, c(panel_list, list(ncol = n_cols)))

final_plot

ggsave("./Figure/Kendrick_mass_defect_plots.pdf",
       plot   = final_plot,
       width  = 1230/96,
       height = 1450/96,
       units  = "in")