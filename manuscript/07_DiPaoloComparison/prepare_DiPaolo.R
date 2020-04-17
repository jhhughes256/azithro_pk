# Prepare digitised azithromycin data from Di Paolo paper
# -----------------------------------------------------------------------------
# This script is designed to prepare digitised plot data to compare with 
# predictions from the Zheng S et al. mrgsolve model.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data digitised from:
#   Di Paolo, A., et al. (2002). "PHARMACOKINETICS OF AZITHROMYCIN IN LUNG 
#   TISSUE, BRONCHIAL WASHING, AND PLASMA IN PATIENTS GIVEN MULTIPLE ORAL DOSES 
#   OF 500 AND 1000 MG DAILY." Pharmacological Research 46(6): 545-550.
# Figure data digitised using open source digitising software Engauge v12.1
# Available at: https://github.com/markummitchell/engauge-digitizer/releases
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Prepare work environment
# Clear objects from R environment
  rm(list = ls(all = TRUE))

# Load packages
  library(tidyverse)
  
# Source external files
# * Read in digitised data
  fig1a_dipaolodig_plasma <- read_csv("DigitisedPlasma_DiPaoloFigure1A.csv", 
    col_types = cols(.default = col_double()))
  fig1b_dipaolodig_plasma <- read_csv("DigitisedPlasma_DiPaoloFigure1B.csv", 
    col_types = cols(.default = col_double()))
  fig1a_dipaolodig_bronchial <- read_csv("DigitisedBronchial_DiPaoloFigure1A.csv", 
    col_types = cols(.default = col_double()))
  fig1b_dipaolodig_bronchial <- read_csv("DigitisedBronchial_DiPaoloFigure1B.csv", 
    col_types = cols(.default = col_double()))
  fig1a_dipaolodig_lung <- read_csv("DigitisedLung_DiPaoloFigure1A.csv", 
    col_types = cols(.default = col_double()))
  fig1b_dipaolodig_lung <- read_csv("DigitisedLung_DiPaoloFigure1B.csv", 
    col_types = cols(.default = col_double()))

# Define nominal times
  nom_times <- c(4, 6, 12, 24, 48, 60, 72, 96, 108, 120, 156, 204)
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare digitised data to be compared
# Digitised data (units: ug/mL)
# * For each figure and sample type
#   + if Figure1A and Plasma data, then generate SD column
#   + Rename x to be Time
#   + Reshape all columns except Time into two columns: column name and value
#   + Convert units from ug/mL to ng/mL
#   + Split column names into Sample (Plasma, Lung etc.) and Metric (MEAN etc.)
#   + Transform each unique Metric into its own column
#   + Round times to closest nominal time
#   + Create SD column and Dose column
#   + Remove processing columns
#   + Calculate geometric mean value for each Sample, Time combination
  fig_names <- c("fig1a", "fig1b")
  sample_names <- c("plasma", "bronchial", "lung")
  figdf <- map2_dfr(rep(fig_names, times = 3), rep(sample_names, each = 2), ~ .x %>%
    paste("dipaolodig", .y, sep = "_") %>%
    get() %>%
    purrr::when(ncol(.) != 3 
      ~ mutate(., Plasma_MeanSD = Plasma_Mean),
      ~ .) %>%
    rename(Time = x) %>%
    pivot_longer(cols = matches("Mean"), names_to = "Sample_Metric", values_to = "value") %>%
    mutate(value = value*1000) %>%
    separate(Sample_Metric, c("Sample", "Metric"), sep = "_") %>%
    pivot_wider(id_cols = c("Sample", "Time"), names_from = "Metric", values_from = "value") %>%
    mutate(whichNom = map_dbl(Time, ~ .x %>%
      magrittr::subtract(nom_times) %>%
      abs() %>%
      {which(. == min(.))})) %>%
    mutate(Time = nom_times[whichNom], SD = abs(MeanSD - Mean)) %>%
    mutate(Dose = case_when(
      .x == "fig1a" ~ 500,
      .x == "fig1b" ~ 1000
    )) %>%
    select(-whichNom, -MeanSD) %>%
    group_by(Sample, Time) %>%
    summarise_all(~ exp(mean(log(.x))))
  )

# Save digitised data
  write_csv(figdf, "DiPaoloFigure1.csv")
  