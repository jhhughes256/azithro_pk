# Prepare digitised azithromycin data from Danesi et al.
# -----------------------------------------------------------------------------
# This script is designed to prepare digitised plot data to compare with 
# predictions from the Zheng S et al. mrgsolve model.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data digitised from:
#   Danesi, R., et al. (2003). "Comparative distribution of azithromycin in 
#   lung tissue of patients given oral daily doses of 500 and 1000 mg." 
#   Journal of Antimicrobial Chemotherapy 51(4): 939-945.
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
  lung_fig1a_danesidig <- read_csv("DigitisedLung_DanesiFigure1A.csv", 
    col_types = cols(.default = col_double()))
  plasma_fig1a_danesidig <- read_csv("DigitisedPlasma_DanesiFigure1A.csv", 
    col_types = cols(.default = col_double()))
  lung_fig1b_danesidig <- read_csv("DigitisedLung_DanesiFigure1B.csv", 
    col_types = cols(.default = col_double()))
  plasma_fig1b_danesidig <- read_csv("DigitisedPlasma_DanesiFigure1B.csv", 
    col_types = cols(.default = col_double()))

# Define digitised & nominal times
  tafd_zero <- 48
  nom_times <- list(
    plas = tafd_zero + c(1, 4, 6, 12, 24, 48, 60, 72, 96, 108, 120, 156, 204),
    lung = tafd_zero + c(6, 12, 60, 108, 156, 204)) 
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare digitised data to be compared
# Digitised data (units: mg/L)
# * For each sample type
#   + Rename x to be Time
#   + Reshape all columns except Time into two columns: Sample and Obs
#   + Convert units from mg/L to ng/mL
#   + Round times to closest nominal time
#   + Remove processing columns
  figdf <- map_dfr(c("lung_fig1a", "plasma_fig1a", "lung_fig1b", "plasma_fig1b"), ~ .x %>%
    paste0("_danesidig") %>%
    get() %>%
    mutate(Time = tafd_zero + x) %>%
    select(-x) %>%
    pivot_longer(cols = -c("Time"), names_to = "Sample_Metric", values_to = "value") %>%
    mutate(value = value*1000, UID = rep(1:(nrow(.)/2), each = 2)) %>%
    separate(Sample_Metric, c("Sample", "Metric"), sep = "_") %>%
    pivot_wider(id_cols = c("Sample", "Time", "UID"), names_from = "Metric", values_from = "value") %>%
    mutate(whichNom = map_dbl(Time, function(.y) .y %>%
      magrittr::subtract(nom_times[[substr(.x, 1, 4)]]) %>%
      abs() %>%
      {which(. == min(.))})) %>%
    mutate(SD = MeanSD - Mean, Time = nom_times[[substr(.x, 1, 4)]][whichNom]) %>%
    select(-whichNom, -UID, -MeanSD) %>%
    group_by(Sample, Time) %>%
    summarise_all(~ exp(mean(log(.x)))) %>%
    mutate(Dose = if_else(str_detect(.x, "fig1a$"), 500, 1000))
  )

# Save digitised data
  write_csv(figdf, "DanesiFigure1.csv")
  