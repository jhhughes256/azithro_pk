# Prepare digitised data from Zheng S et al.
# ------------------------------------------------------------------------------
# This script aims to prepare digitised data to compare with the model from
# the same manuscript.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data sourced from:
#   Songmao Zheng, Peter Matzneller, Markus Zeitlinger, Stephan Schmidt, 
#   Development of a Population Pharmacokinetic Model Characterizing the Tissue
#   Distribution of Azithromycin in Healthy Subjects, Antimicrobial Agents and 
#   Chemotherapy Oct 2014, 58 (11) 6675-6684; DOI: 10.1128/AAC.02904-14
# Data digitised using open source digitising software Engauge v12.1
# Available at: https://github.com/markummitchell/engauge-digitizer/releases
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare work environment
# Clear objects from R environment
  rm(list = ls(all = TRUE))

# Load packages
  library(tidyverse)

# Read in digitised data
  plas_zhengdig <- read_csv("digitisedObs_PLAS_ZhengPlot.csv", 
    col_types = cols(.default = col_double())) %>%
    rename(Time = x, PLAS_MEAN = PLAS_OBS)
  musc_zhengdig <- read_csv("digitisedObs_MUSC_ZhengPlot.csv", 
    col_types = cols(.default = col_double())) %>%
      rename(Time = x, MUSC_MEAN = MUSC_OBS)
  subc_zhengdig <- read_csv("digitisedObs_SUBC_ZhengPlot.csv", 
    col_types = cols(.default = col_double())) %>%
      rename(Time = x, SUBC_MEAN = SUBC_OBS)
  pmlc_zhengdig <- read_csv("digitisedObs_PML_ZhengPlot.csv", 
    col_types = cols(.default = col_double())) %>%
      rename(Time = x, PMLC_MEAN = PMLC_OBS)
  
# Define nominal times
  nom_times <- list(
    PLAS = c(
      0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 6, 8,  # day 1
      (48 + c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 6, 8)),  # day 3
      (96 + c(0, 2, 4)), # day 5
      (216 + c(0, 2, 4))), # day 10
    MUSC = c(
      0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5,  # day 1
      (48 + c(0, 0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5)),  # day 3
      (96 + c(1, 3)), # day 5
      (216 + c(1, 3))), # day 10
    PMLC = c(
      0, 2, 6, 10,  # day 1
      (48 + c(0, 2, 6, 10)),  # day 3
      96, # day 5
      216)) # day 10
  nom_times$SUBC <- nom_times$MUSC
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare digitised data for plotting
# Digitised data (units: ng/mL)
# * Bind data together
# * Reshape all columns except Sample into two columns: column name and value
# * Split column names into Sample (SERUM, AM, etc.) and Metric (MEAN etc.)
# * Transform each unique Metric into its own column
# * Create SE column and remove MEAN+SE column (originally named SE)
# * For each sample type
#   + Determine which nominal time is closest to the digitised time
#   + Change Time to nominal time for that Sample
# * For each sample type and nominal time
#   + Find the geometric mean of the MEAN and SD values (Geometric mean is 
#     used as data was digitised on a log-scale, i.e. clicking a pixel above
#     the point leads to greater error than a pixel below)
  figdf <- bind_rows(plas_zhengdig, musc_zhengdig, subc_zhengdig, pmlc_zhengdig) %>%
    pivot_longer(cols = matches("PLAS|MUSC|SUBC|PMLC"), names_to = "Sample_Metric") %>%
    separate(Sample_Metric, c("Sample", "Metric"), sep = "_") %>%
    filter(!is.na(value)) %>%
    pivot_wider(id_cols = c("Sample", "Time"), names_from = "Metric") %>%
    mutate(SE = abs(SE - MEAN)) %>%
    group_by(Sample) %>%
    mutate(whichNom = map_dbl(Time, ~ which( abs(.x - nom_times[[unique(Sample)]]) == min(abs(.x - nom_times[[unique(Sample)]])) ))) %>%
    mutate(Time = nom_times[[unique(Sample)]][whichNom]) %>%
    select(-whichNom) %>%
    group_by(Sample, Time) %>%
    summarise_all(~ exp(mean(log(.x))))
  
# Write prepared data to .csv file
  write_csv(figdf, "ZhengFigure3.csv")
  