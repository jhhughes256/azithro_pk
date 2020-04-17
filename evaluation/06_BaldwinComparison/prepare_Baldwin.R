# Prepare digitised azithromycin data from Baldwin paper
# -----------------------------------------------------------------------------
# This script is designed to prepare digitised plot data to compare with 
# predictions from the Zheng S et al. mrgsolve model.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data digitised from:
#   Baldwin, D. R., et al. (1990). "Azithromycin concentrations at the sites of 
#   pulmonary infection." Eur Respir J 3(8): 886-890.
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
  fig2_baldwindig <- read_csv("Digitised_BaldwinFigure2.csv", 
    col_types = cols(.default = col_double()))

# Define digitised & nominal times
  dig_times <- c(1, 2, 3, 4, 5)
  nom_times <- c(12, 24, 48, 72, 96)
  
# Define number of patients from whom samples were taken at each nominal time
  n_vals <- c(4, 4, 4, 6, 4)
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare digitised data to be compared
# Digitised data (units: mg/L)
# * Rename x to Time and convert values to ng/mL
# * Reshape all columns except Time into two columns: column name and value
# * Split column names into Form (IR, ER, etc.) and Metric (MEAN etc.)
# * Transform each unique Metric into its own column
# * Remove imputed times created by digitising software
#   + This is done by determining the digitised time for each set of columns
#   + Subtract this from actual digitised time to get difference
#   + Digitised values for each sample are: Serum < 0, 0 < ELF < 0.2, AM > 0.2
# * Create SD column by calculating SEM then multiply by sqrt(N)
# * Remove processing columns
# * Calculate mean value for each Sample, Time combination
  figdf <- fig2_baldwindig %>%
    rename(Time = x) %>%
    mutate_at(vars(matches("MEAN")), ~ .x * 1000) %>%
    pivot_longer(cols = matches("MEAN"), names_to = "Sample_Metric", values_to = "value") %>%
    separate(Sample_Metric, c("Sample", "Metric"), sep = "_") %>%
    pivot_wider(id_cols = c("Sample", "Time"), names_from = "Metric", values_from = "value") %>%
    mutate(digTime = round(Time, 0), deltaDigTime = Time - digTime) %>%
    filter((Sample == "SERUM" & deltaDigTime < 0) | 
      (Sample == "ELF" & deltaDigTime > 0 & deltaDigTime < 0.2) | 
      (Sample == "AM" & deltaDigTime > 0.2)) %>%
    mutate(Time = nom_times[digTime]) %>%
    mutate(SD = abs(MEANPLUSSEM - MEAN)*sqrt(n_vals[digTime])) %>%
    select(-MEANPLUSSEM, -digTime, -deltaDigTime) %>%
    group_by(Sample, Time) %>%
    summarise_all(mean)

# Save digitised data
  write_csv(figdf, "BaldwinFigure2.csv")
  