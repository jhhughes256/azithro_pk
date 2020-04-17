# Prepare digitised azithromycin data from A0661112 Study Data
# -----------------------------------------------------------------------------
# This script is designed to prepare digitised plot data to compare with 
# predictions from the Zheng S et al. mrgsolve model.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data digitised from:
#   Clinical Study A0661112
#   A Phase 1, Open-Label, Randomized, Parallel Group Study to Compare Serum
#   and White Blood Cell Pharmacokinetics of a 2 Gram Single Dose Sustained 
#   Release Azithromycin Formulation and a 3-Day Regimen of Immediate Release 
#   Azithromycin Commercial Tablet Formulation in Healthy Adult Subjects 
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
  fig1_studydig <- read_csv("Digitised_A0661112FigureS1.csv", 
    col_types = cols(.default = col_double()))
  fig2_studydig <- read_csv("Digitised_A0661112FigureS2.csv", 
    col_types = cols(.default = col_double()))
  fig3_studydig <- read_csv("Digitised_A0661112FigureS3.csv", 
    col_types = cols(.default = col_double()))

# Define digitised & nominal times
  nom_times <- list(
    fig1 = c(
      0.5, 1, 2, 3, 4, 6, 8, 12, 24, 48,  # day 1, 2
      (48 + c(0.5, 1, 2, 3, 4, 6, 8, 12)),  # day 3
      72, 96, 120),  # day 4, 5, 6
    fig2 = c(2, 4, 8, 12, 24, 28, 36, 48, 52, 60, 72, 96, 120))
  nom_times$fig3 <- nom_times$fig2
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare digitised data to be compared
# Digitised data (units: ug/mL)
# * For each figure and sample type
#   + Rename x to be Time
#   + Reshape all columns except Time into two columns: Form and Mean
#   + Filter out concentrations below zero
#   + Convert units from ug/mL to ng/mL
#   + Round times to closest nominal time
#   + Create SD column and Sample column
#   + Remove processing columns
#   + Calculate mean value for each Form, Sample, Time combination
  figdf <- map_dfr(c("fig1", "fig2", "fig3"), ~ .x %>%
    paste0("_studydig") %>%
    get() %>%
    rename(Time = x) %>%
    pivot_longer(cols = c("SR", "IR"), names_to = "Form", values_to = "Mean") %>%
    filter(Mean > 0) %>%
    mutate(Mean = Mean*1000) %>%
    mutate(whichNom = map_dbl(Time, function(.y) .y %>%
      magrittr::subtract(nom_times[[.x]]) %>%
      abs() %>%
      {which(. == min(.))})) %>%
    mutate(Time = nom_times[[.x]][whichNom]) %>%
    mutate(Sample = case_when(
      .x == "fig1" ~ "Serum",
      .x == "fig2" ~ "MNL",
      .x == "fig3" ~ "PML"
    )) %>%
    select(-whichNom) %>%
    group_by(Form, Sample, Time) %>%
    summarise_all(mean)
  )

# Save digitised data
  write_csv(figdf, "A0661112FigureS123.csv")
  