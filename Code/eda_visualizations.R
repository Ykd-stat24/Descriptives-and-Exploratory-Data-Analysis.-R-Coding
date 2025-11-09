# Exploratory visualizations for LD 2025 dataset
# ------------------------------------------------
# This script loads the processed pairwise discrimination dataset and
# generates a suite of exploratory figures required for QC and reporting.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(janitor)
  library(fs)
  library(survival)
  library(survminer)
  library(ggridges)
})

# ---------------------------------------------------------------------
# 1. Load and clean data ------------------------------------------------
# ---------------------------------------------------------------------

raw_path <- "Processed Data/W LD Data 11.08 PD_processed.csv"
if (!file.exists(raw_path)) {
  stop("Data file not found at ", raw_path)
}

raw_data <- readr::read_csv(raw_path, na = c("", "NA", "NaN"), show_col_types = FALSE)

clean_data <- raw_data %>% janitor::clean_names()

required_columns <- c(
  "schedule_name", "schedule_run_date", "animal_id", "day",
  "genotype", "end_summary_percent_correct_1",
  "end_summary_trials_completed_1", "correct_touch_latency_mean",
  "incorrect_touch_latency_mean"
)

missing_columns <- setdiff(required_columns, names(clean_data))
if (length(missing_columns) > 0) {
  stop(
    "Missing expected columns after cleaning: ",
    paste(missing_columns, collapse = ", ")
  )
}

ld_data <- clean_data %>%
  mutate(
    schedule_name = as.character(schedule_name),
    animal_id = str_to_lower(animal_id),
    cohort = str_to_upper(str_extract(animal_id, "^c\\d+")),
    cohort = forcats::fct_inorder(cohort),
    genotype = forcats::fct_relevel(as.factor(genotype)),
    day = as.integer(day),
    session_date = parse_date_time(schedule_run_date, orders = c("mdy HM", "mdy HMS"), quiet = TRUE)
  )

ld_data <- ld_data %>%
  mutate(
    phase = case_when(
      str_detect(schedule_name, "Initial Touch") ~ "Initial Touch",
      str_detect(schedule_name, "Must Touch") & str_detect(schedule_name, "Pairwise") ~ "Pairwise Must Touch",
      str_detect(schedule_name, "Must Touch") ~ "Must Touch",
      str_detect(schedule_name, "Punish Incorrect") & str_detect(schedule_name, "Pairwise") ~ "Pairwise Punish Incorrect",
      str_detect(schedule_name, "Punish Incorrect") ~ "Punish Incorrect",
      str_detect(schedule_name, "1 choice reversal") ~ "1 Choice Reversal",
      str_detect(schedule_name, "1 choice") ~ "1 Choice",
      str_detect(schedule_name, "Pairwise Discrimination") & str_detect(schedule_name, "Reversal") ~ "Pairwise Discrimination Reversal",
      str_detect(schedule_name, "Pairwise Discrimination") ~ "Pairwise Discrimination",
      TRUE ~ NA_character_
    ),
    phase = factor(phase, levels = c(
      "Initial Touch",
      "Must Touch",
      "Punish Incorrect",
      "1 Choice",
      "1 Choice Reversal",
      "Pairwise Must Touch",
      "Pairwise Punish Incorrect",
      "Pairwise Discrimination",
      "Pairwise Discrimination Reversal"
    )),
      percent_correct = end_summary_percent_correct_1,
      trials_completed = end_summary_trials_completed_1,
    correct_latency = correct_touch_latency_mean,
    incorrect_latency = incorrect_touch_latency_mean
  )

# Drop sessions that cannot be classified into a phase
ld_data <- ld_data %>% filter(!is.na(phase))

# Default colour palette for genotypes ---------------------------------

genotype_palette <- c(
  "+" = "#1b9e77",
  "+/-" = "#d95f02",
  "-" = "#7570b3"
)

genotype_palette <- genotype_palette[names(genotype_palette) %in% levels(ld_data$genotype)]

# Ensure figures directory exists
fs::dir_create("figures")

# ---------------------------------------------------------------------
# 2. Percent correct trajectories --------------------------------------
# ---------------------------------------------------------------------

percent_correct_plot <- ld_data %>%
  ggplot(aes(x = day, y = percent_correct, colour = genotype, group = interaction(animal_id, genotype))) +
  geom_line(alpha = 0.6) +
  geom_point(alpha = 0.6, size = 1) +
  facet_wrap(~ cohort, ncol = 2, scales = "free_x") +
  scale_colour_manual(values = genotype_palette, na.translate = FALSE) +
  labs(
    title = "Percent correct across training days",
    x = "Session day",
    y = "Percent correct",
    colour = "Genotype"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.spacing = grid::unit(1, "lines"))

ggsave(
  filename = "figures/eda_percent_correct_by_cohort.png",
  plot = percent_correct_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# ---------------------------------------------------------------------
# 3. Trials-to-criterion and Kaplan-Meier curves -----------------------
# ---------------------------------------------------------------------

criterion_threshold <- 80

criterion_summary <- ld_data %>%
  arrange(animal_id, phase, session_date, day) %>%
  group_by(animal_id, phase) %>%
  mutate(
    session_in_phase = row_number(),
    attained = percent_correct >= criterion_threshold
  ) %>%
  summarise(
    genotype = dplyr::first(genotype),
    cohort = dplyr::first(cohort),
    sessions_observed = dplyr::n(),
    reached_criterion = any(attained, na.rm = TRUE),
    attainment_session = {
      idx <- which(attained)
      if (length(idx) > 0) {
        min(session_in_phase[idx], na.rm = TRUE)
      } else {
        sessions_observed
      }
    },
    .groups = "drop"
  ) %>%
  mutate(
    attainment_session = replace_na(attainment_session, sessions_observed)
  )

surv_object <- Surv(time = criterion_summary$attainment_session, event = criterion_summary$reached_criterion)

km_fit <- survfit(surv_object ~ genotype + phase, data = criterion_summary)

km_plot <- ggsurvplot_facet(
  km_fit,
  data = criterion_summary,
  facet.by = "phase",
  fun = "event",
  legend.title = "Genotype",
  ggtheme = theme_minimal(base_size = 11),
  palette = genotype_palette,
  break.time.by = 2,
  xlab = "Sessions within phase",
  ylab = "Cumulative proportion reaching criterion"
)

km_plot$plot <- km_plot$plot +
  theme(strip.text = element_text(face = "bold"))

ggsave(
  filename = "figures/eda_trials_to_criterion_km.png",
  plot = km_plot$plot,
  width = 10,
  height = 8,
  dpi = 300
)

# ---------------------------------------------------------------------
# 4. Latency distributions (ridge plots) -------------------------------
# ---------------------------------------------------------------------

latency_long <- ld_data %>%
  select(animal_id, genotype, cohort, phase, day, correct_latency, incorrect_latency) %>%
  pivot_longer(
    cols = c(correct_latency, incorrect_latency),
    names_to = "latency_type",
    values_to = "latency"
  ) %>%
  mutate(
    latency_type = recode(latency_type,
                          correct_latency = "Correct",
                          incorrect_latency = "Blank/Incorrect"),
    latency_type = factor(latency_type, levels = c("Correct", "Blank/Incorrect"))
  ) %>%
  drop_na(latency)

latency_ridge_plot <- latency_long %>%
  ggplot(aes(x = latency, y = phase, fill = genotype)) +
  geom_density_ridges(alpha = 0.7, colour = NA, scale = 1.1) +
  facet_wrap(~ latency_type, scales = "free_x") +
  scale_fill_manual(values = genotype_palette, na.translate = FALSE) +
  labs(
    title = "Latency distributions by phase and genotype",
    x = "Latency (s)",
    y = "Phase",
    fill = "Genotype"
  ) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

ggsave(
  filename = "figures/eda_latency_ridges.png",
  plot = latency_ridge_plot,
  width = 10,
  height = 7,
  dpi = 300
)

# ---------------------------------------------------------------------
# 5. Session-to-session transition heatmaps ----------------------------
# ---------------------------------------------------------------------

heatmap_plot <- ld_data %>%
  mutate(
    phase = fct_drop(phase),
    animal_label = str_to_upper(animal_id)
  ) %>%
  ggplot(aes(x = day, y = fct_rev(phase), fill = percent_correct)) +
  geom_tile(colour = "white", linewidth = 0.1) +
  facet_wrap(~ animal_label, ncol = 4) +
  scale_fill_viridis_c(option = "C", limits = c(0, 100), na.value = "grey90") +
  labs(
    title = "Session-to-session percent correct heatmaps",
    x = "Session day",
    y = "Phase",
    fill = "% correct"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank()
  )

ggsave(
  filename = "figures/eda_percent_correct_heatmaps.png",
  plot = heatmap_plot,
  width = 12,
  height = 9,
  dpi = 300
)

# ---------------------------------------------------------------------
# 6. Trials-to-criterion summary table ---------------------------------
# ---------------------------------------------------------------------

criterion_summary_table <- criterion_summary %>%
  select(animal_id, cohort, phase, genotype, attainment_session, reached_criterion)

readr::write_csv(criterion_summary_table, "figures/eda_trials_to_criterion_summary.csv")

message("EDA figures exported to the figures/ directory.")
