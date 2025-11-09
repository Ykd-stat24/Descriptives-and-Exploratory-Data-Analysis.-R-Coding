#' Location Discrimination baseline mixed-effects modelling
#'
#' This script prepares session-level outcomes for the LD, PD, and PDR phases
#' and fits baseline mixed-effects models for accuracy and latency outcomes.
#' It expects the processed CSV exports located in the `Processed Data/` folder.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(glmmTMB)
  library(lme4)
  library(DHARMa)
  library(performance)
  library(broom.mixed)
  library(ggeffects)
})

# -----------------------------------------------------------------------------
# Helper functions -------------------------------------------------------------
# -----------------------------------------------------------------------------

read_ld_data <- function(path, cohort_label) {
  raw <- readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(
      .row_id = dplyr::row_number(),
      Cohort = cohort_label,
      Phase = dplyr::case_when(
        stringr::str_detect(`Schedule name`, regex("LD 1 choice", ignore_case = TRUE)) ~ "LD Choice",
        stringr::str_detect(`Schedule name`, regex("Pairwise Discrimination v3 - Reversal", ignore_case = TRUE)) ~ "PDR",
        stringr::str_detect(`Schedule name`, regex("Pairwise Discrimination", ignore_case = TRUE)) ~ "PD",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(Phase))

  trial_correct_cols <- tidyselect::matches("^Trial Analysis - No\\. Correct \\\(")
  correct_latency_cols <- tidyselect::matches("^Trial Analysis - Correct Image Response Latency \\\(")
  reward_latency_cols <- tidyselect::matches("^Trial Analysis - Reward Collection Latency \\\(")

  correct_counts <- raw %>%
    dplyr::select(.row_id, !!trial_correct_cols) %>%
    dplyr::mutate(dplyr::across(-.row_id, readr::parse_number)) %>%
    dplyr::mutate(
      n_correct = rowSums(dplyr::across(-.row_id), na.rm = TRUE),
      n_trials = rowSums(!is.na(dplyr::across(-.row_id)))
    ) %>%
    dplyr::mutate(n_incorrect = pmax(n_trials - n_correct, 0)) %>%
    dplyr::select(.row_id, n_correct, n_incorrect, n_trials)

  correct_latency <- raw %>%
    dplyr::select(.row_id, !!correct_latency_cols) %>%
    dplyr::mutate(dplyr::across(-.row_id, readr::parse_number)) %>%
    dplyr::mutate(correct_latency_mean = rowMeans(dplyr::across(-.row_id), na.rm = TRUE)) %>%
    dplyr::mutate(correct_latency_mean = dplyr::na_if(correct_latency_mean, Inf)) %>%
    dplyr::select(.row_id, correct_latency_mean)

  reward_latency <- raw %>%
    dplyr::select(.row_id, !!reward_latency_cols) %>%
    dplyr::mutate(dplyr::across(-.row_id, readr::parse_number)) %>%
    dplyr::mutate(reward_latency_mean = rowMeans(dplyr::across(-.row_id), na.rm = TRUE)) %>%
    dplyr::mutate(reward_latency_mean = dplyr::na_if(reward_latency_mean, Inf)) %>%
    dplyr::select(.row_id, reward_latency_mean)

  cleaned <- raw %>%
    dplyr::mutate(
      schedule_datetime = lubridate::parse_date_time(
        `Schedule run date`,
        orders = c("mdY HMS", "mdY HM", "m/d/Y H:M")
      ),
      trials_completed = readr::parse_number(`End Summary - Trials Completed (1)`),
      percent_correct = readr::parse_number(`End Summary - Percentage Correct (1)`),
      AnimalID = `Animal ID`
    ) %>%
    dplyr::left_join(correct_counts, by = ".row_id") %>%
    dplyr::left_join(correct_latency, by = ".row_id") %>%
    dplyr::left_join(reward_latency, by = ".row_id") %>%
    dplyr::mutate(
      n_trials = dplyr::coalesce(n_trials, trials_completed),
      accuracy = dplyr::case_when(
        !is.na(n_trials) & n_trials > 0 ~ n_correct / n_trials,
        !is.na(percent_correct) ~ percent_correct / 100,
        TRUE ~ NA_real_
      ),
      accuracy = pmin(pmax(accuracy, 1e-4), 1 - 1e-4),
      logit_accuracy = log(accuracy / (1 - accuracy)),
      Genotype = dplyr::case_when(
        Genotype %in% c("+", "+/+", "+/+") ~ "B6",
        Genotype %in% c("-", "-/ -", "-/-") ~ "C3H",
        Genotype %in% c("+/-", "-/+") ~ "Het",
        TRUE ~ Genotype
      ),
      Genotype = factor(Genotype),
      Cohort = factor(Cohort)
    ) %>%
    dplyr::arrange(AnimalID, Phase, schedule_datetime, readr::parse_number(Day)) %>%
    dplyr::group_by(AnimalID, Phase) %>%
    dplyr::mutate(
      session_index = dplyr::row_number(),
      session_scaled = as.numeric(scale(session_index))
    ) %>%
    dplyr::ungroup()

  cleaned
}

# -----------------------------------------------------------------------------
# Load and combine the Winter/Fall processed exports --------------------------
# -----------------------------------------------------------------------------

ld_winter <- read_ld_data("Processed Data/W LD Data 11.08 All_processed.csv", "Winter")
ld_fall   <- read_ld_data("Processed Data/F LD Data 11.08 All_processed.csv", "Fall")

ld_sessions <- dplyr::bind_rows(ld_winter, ld_fall) %>%
  dplyr::filter(!is.na(accuracy), !is.na(Genotype), !is.na(Cohort)) %>%
  dplyr::mutate(
    Phase = factor(Phase, levels = c("LD Choice", "PD", "PDR")),
    session_scaled = dplyr::if_else(is.na(session_scaled), 0, session_scaled)
  )

readr::write_csv(ld_sessions, "Processed Data/ld_sessions_model_input.csv")

# -----------------------------------------------------------------------------
# Mixed-effects models --------------------------------------------------------
# -----------------------------------------------------------------------------

# Accuracy (binomial mixed model)
accuracy_formula_add <- cbind(n_correct, n_incorrect) ~ Genotype + Cohort + session_scaled + Phase + (1 + session_scaled | AnimalID)
accuracy_formula_full <- cbind(n_correct, n_incorrect) ~ Genotype * Cohort * session_scaled + Phase + (1 + session_scaled | AnimalID)

accuracy_model_add <- glmmTMB::glmmTMB(
  formula = accuracy_formula_add,
  family  = binomial,
  data    = ld_sessions
)

accuracy_model_full <- glmmTMB::glmmTMB(
  formula = accuracy_formula_full,
  family  = binomial,
  data    = ld_sessions
)

accuracy_lrtest <- anova(accuracy_model_add, accuracy_model_full)
accuracy_r2     <- performance::r2(accuracy_model_full)

# Latency (gamma mixed model on correct response latency)
latency_data <- ld_sessions %>%
  dplyr::filter(!is.na(correct_latency_mean), correct_latency_mean > 0)

latency_formula_add <- correct_latency_mean ~ Genotype + Cohort + session_scaled + Phase + (1 + session_scaled | AnimalID)
latency_formula_full <- correct_latency_mean ~ Genotype * Cohort * session_scaled + Phase + (1 + session_scaled | AnimalID)

latency_model_add <- glmmTMB::glmmTMB(
  formula = latency_formula_add,
  family  = Gamma(link = "log"),
  data    = latency_data
)

latency_model_full <- glmmTMB::glmmTMB(
  formula = latency_formula_full,
  family  = Gamma(link = "log"),
  data    = latency_data
)

latency_lrtest <- anova(latency_model_add, latency_model_full)
latency_r2     <- performance::r2(latency_model_full)

# -----------------------------------------------------------------------------
# Diagnostic checks -----------------------------------------------------------
# -----------------------------------------------------------------------------

accuracy_res <- DHARMa::simulateResiduals(accuracy_model_full)
plot(accuracy_res)
DHARMa::testDispersion(accuracy_res)
DHARMa::testZeroInflation(accuracy_res)

latency_res <- DHARMa::simulateResiduals(latency_model_full)
plot(latency_res)
DHARMa::testDispersion(latency_res)

# -----------------------------------------------------------------------------
# Summaries & effect plots ----------------------------------------------------
# -----------------------------------------------------------------------------

accuracy_tidy <- broom.mixed::tidy(accuracy_model_full, effects = "fixed", conf.int = TRUE, exponentiate = TRUE)
latency_tidy  <- broom.mixed::tidy(latency_model_full, effects = "fixed", conf.int = TRUE)

readr::write_csv(accuracy_tidy, "Processed Data/accuracy_model_fixed_effects.csv")
readr::write_csv(latency_tidy,  "Processed Data/latency_model_fixed_effects.csv")

accuracy_effects <- ggeffects::ggpredict(
  accuracy_model_full,
  terms = c("session_scaled [all]", "Genotype", "Cohort"),
  type  = "re" # conditional on random effects = population level
)

latency_effects <- ggeffects::ggpredict(
  latency_model_full,
  terms = c("session_scaled [all]", "Genotype", "Cohort"),
  type  = "re"
)

acc_plot <- plot(accuracy_effects) +
  ggplot2::labs(
    x = "Session (scaled)",
    y = "Predicted accuracy (probability of correct)",
    title = "Predicted accuracy across sessions"
  ) +
  ggplot2::theme_minimal()

ggplot2::ggsave("Processed Data/accuracy_effects_plot.png", acc_plot, width = 8, height = 5, dpi = 300)

lat_plot <- plot(latency_effects) +
  ggplot2::labs(
    x = "Session (scaled)",
    y = "Predicted correct response latency (s)",
    title = "Predicted latency across sessions"
  ) +
  ggplot2::theme_minimal()

ggplot2::ggsave("Processed Data/latency_effects_plot.png", lat_plot, width = 8, height = 5, dpi = 300)

# Display key outputs in the console when the script is sourced ----------------

print(dplyr::glimpse(ld_sessions))
print(accuracy_lrtest)
print(accuracy_r2)
print(latency_lrtest)
print(latency_r2)
print(accuracy_tidy)
print(latency_tidy)
