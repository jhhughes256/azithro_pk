# Azithromycin Population Pharmacokinetic Model
# ------------------------------------------------------------------------------
# This script is designed to generate a population for simulation of azithromycin 
# concentrations using the regimen described in recent literature on COVID-19:
#   Gautret et al. Hydroxychloroquine and azithromycin as a treatment of 
#   COVID-19: results of an open-label non-randomized clinical trial, 
#   International Journal of Antimicrobial Agents, Mar 2020; 
#   DOI: 10.1016/j/ijantimicag.2020.105949
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load package libraries
  library(tidyverse)

# Define population
  nid <- 1000
  dose_amt <- c(500, rep(250, times = 4))
  dose_times <- seq(0, 96, by = 24)
  samp_times <- seq(0, 120, by = 0.5)
  
# Create population dataset
  dosedf <- data.frame(
    C = ".",
    ID = rep(1:nid, each = length(dose_times)),
    TIME = rep(dose_times, times = nid),
    AMT = rep(dose_amt, times = nid),
    CMT = 1,
    EVID = 1,
    MDV = 1
  )
  sampdf <- data.frame(
    C = ".",
    ID = rep(1:nid, each = length(samp_times)),
    TIME = rep(samp_times, times = nid),
    AMT = 0,
    CMT = NA,
    EVID = 0,
    MDV = 0
  )
  popdf <- rbind(dosedf,
      mutate(sampdf, CMT = 2),
      mutate(sampdf, CMT = 3),
      mutate(sampdf, CMT = 4),
      mutate(sampdf, CMT = 9)) %>%
    arrange(ID, TIME, CMT) %>%
    mutate(DV = 0, DVID = CMT)

# Save to file
  write_csv(popdf, "evaluation/02_NONMEMValidation/AZPop_30MAR2020.csv")
  
  
