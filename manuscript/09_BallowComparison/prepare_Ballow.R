# Prepare digitised azithromycin data from Ballow et al.
# -----------------------------------------------------------------------------
# This script is designed to prepare digitised plot data to compare with 
# predictions from the Zheng S et al. mrgsolve model.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data digitised from:
#   Ballow, C. H., et al. (1998). "Pharmacokinetics of Oral Azithromycin in 
#   Serum, Urine, Polymorphonuclear Leucocytes and Inflammatory vs 
#   Non-Inflammatory Skin Blisters in Healthy Volunteers." Clinical Drug 
#   Investigation 15(2): 159-167.
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
  pml_ballowdig <- read_csv("DigitisedPMN_BallowFigure1.csv", 
    col_types = cols(.default = col_double()))
  serum_ballowdig <- read_csv("DigitisedSerum_BallowFigure1.csv", 
    col_types = cols(.default = col_double()))

# Define digitised & nominal times
  nom_times <- list(
    pml = c(
      4, (24 + c(0, 4)), (48 + c(0, 8)), (72 + c(0, 8)), 96, # day 1, 2, 3, 4
      120, 144, 168, 192, 216, 240, 264),  # day 5, 6, 7, 8, 9, 10, 11
    serum = c(
      0.25, 0.5, 0.75, 1.0, 1.5, 2, 3, 4, 5, 6, 8, 12, 16, 24,  # day 1
      (24 + c(2, 12)), (48 + c(2, 12)), (72 + c(2, 12)),  # day 2, 3, 4
      (96 + c(0.25, 0.5, 0.75, 1.0, 1.5, 2, 3, 4, 5, 6, 8, 12, 16, 24)),  # day 5
      144, 168, 192, 216))  # day 6, 7, 8, 9
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare digitised data to be compared
# Digitised data (units: mg/L)
# * For each sample type
#   + Rename x to be Time
#   + Reshape all columns except Time into two columns: Sample and Obs
#   + Convert units from mg/L to ng/mL
#   + Round times to closest nominal time
#   + Remove processing columns
  figdf <- map_dfr(c("pml", "serum"), ~ .x %>%
    paste0("_ballowdig") %>%
    get() %>%
    rename(Time = x) %>%
    pivot_longer(cols = -c("Time"), names_to = "Sample", values_to = "Obs") %>%
    mutate(Obs = Obs*1000) %>%
    mutate(whichNom = map_dbl(Time, function(.y) .y %>%
      magrittr::subtract(nom_times[[.x]]) %>%
      abs() %>%
      {which(. == min(.))})) %>%
    mutate(Time = nom_times[[.x]][whichNom]) %>%
    select(-whichNom)
  )

# Save digitised data
  write_csv(figdf, "BallowFigure1.csv")
  