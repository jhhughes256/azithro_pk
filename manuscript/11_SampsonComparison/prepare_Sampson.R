# Prepare digitised azithromycin data from Sampson et al.
# -----------------------------------------------------------------------------
# This script is designed to prepare digitised plot data to compare with 
# predictions from the Zheng S et al. mrgsolve model.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data digitised from:
#   Sampson, M. R., et al. (2014). "Population Pharmacokinetics of Azithromycin 
#   in Whole Blood, Peripheral Blood Mononuclear Cells, and Polymorphonuclear 
#   Cells in Healthy Adults." CPT: pharmacometrics & systems 
#   pharmacology 3(3): 103.
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
  pml_fig1a_sampsondig <- read_csv("DigitisedPML_SampsonFigure1A.csv", 
    col_types = cols(.default = col_double()))
  mnl_fig1a_sampsondig <- read_csv("DigitisedMNL_SampsonFigure1A.csv", 
    col_types = cols(.default = col_double()))
  pml_fig1b_sampsondig <- read_csv("DigitisedPML_SampsonFigure1B.csv", 
    col_types = cols(.default = col_double()))
  mnl_fig1b_sampsondig <- read_csv("DigitisedMNL_SampsonFigure1B.csv", 
    col_types = cols(.default = col_double()))

# Define nominal times
  nom_times <- c(1, 2, 3, 4, 6, 9, 12, 16, 24, 48, 96, 144, 240, 336, 504)
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare digitised data to be compared
# Digitised data (units: mg/L)
# * For each sample type
#   + Rename x to be Time
#   + Reshape all columns except Time into two columns: Sample and Obs
#   + Convert units from mg/L to ng/mL
#   + Round times to closest nominal time
#   + Remove processing columns
  figdf <- map_dfr(c("pml_fig1a", "mnl_fig1a", "pml_fig1b", "mnl_fig1b"), ~ .x %>%
    paste0("_sampsondig") %>%
    get() %>%
    rename(Time = x) %>%
    pivot_longer(cols = -c("Time"), names_to = "Sample_Dose_Metric", values_to = "value") %>%
    mutate(value = value*1000, UID = rep(1:(nrow(.)/2), each = 2)) %>%
    separate(Sample_Dose_Metric, c("Sample", "Dose", "Metric"), sep = "_") %>%
    pivot_wider(id_cols = c("Sample", "Dose", "Time", "UID"), names_from = "Metric", values_from = "value") %>%
    purrr::when("MeanSD" %in% names(.) ~ rename(., MeanSE = MeanSD), ~ .) %>%
    mutate(whichNom = map_dbl(Time, function(.y) .y %>%
      magrittr::subtract(nom_times) %>%
      abs() %>%
      {which(. == min(.))})) %>%
    mutate(SE = abs(MeanSE - Mean), Time = nom_times[whichNom]) %>%
    select(-whichNom, -UID, -MeanSE) %>%
    group_by(Sample, Dose, Time) %>%
    summarise_all(~ exp(mean(log(.x))))
  ) %>% arrange(Sample, Dose, Time)

# Save digitised data
  write_csv(figdf, "SampsonFigure1AB.csv")
  