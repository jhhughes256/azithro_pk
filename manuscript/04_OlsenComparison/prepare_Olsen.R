# Prepare digitised azithromycin data from Olsen paper
# -----------------------------------------------------------------------------
# This script is designed to compare digitised plot data with table data and 
# produce a combined dataset to compare with predictions from the Zheng S et 
# al. mrgsolve model.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data digitised from:
#   Olsen, K. M., et al. (1996). "Intrapulmonary pharmacokinetics of 
#   azithromycin in healthy volunteers given five oral doses." Antimicrobial 
#   Agents and Chemotherapy 40(11): 2582.
# Figure data digitised using open source digitising software Engauge v12.1
# Available at: https://github.com/markummitchell/engauge-digitizer/releases
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Prepare work environment
# Clear objects from R environment
  rm(list = ls(all = TRUE))

# Load packages
  library(tidyverse)
  
# Source external files
# * Read in table data
  source("tribble_OlsenTable1.R")  # `olsendf`
 
# * Read in digitised data
  fig1_olsendig <- read_csv("Digitised_OlsenFigure1.csv", 
    col_types = cols(.default = col_double()))
  fig2_olsendig <- read_csv("Digitised_OlsenFigure2.csv", 
    col_types = cols(.default = col_double()))
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare digitised data to be compared
# Table data (units: ug/mL)
# * Convert to ng/mL
# * Reshape all columns except Sample into two columns: column name and value
# * Split column names into Metric (MN, SD) and Time
# * Transform each unique Metric into its own column
# * Rename column and turn time into numeric class
  tabdf <- olsendf %>%
    mutate_at(vars(matches("\\d+$")), ~ .x * 1000) %>%
    pivot_longer(cols = matches("\\d+$"), names_to = "Metric_Time") %>%
    separate(Metric_Time, c("Metric", "Time"), sep = "_") %>%
    pivot_wider(id_cols = c("Sample", "Time"), names_from = "Metric") %>%
    rename(MEAN = MN) %>%
    mutate(Time = as.double(Time))
  
# Digitised data (units: ug/mL)
# * For each dataset
#   + Rename x to Time and convert values to ng/mL
#   + Reshape all columns except Sample into two columns: column name and value
#   + Split column names into Sample (SERUM, AM, etc.) and Metric (MEAN etc.)
#   + Transform each unique Metric into its own column
#   + Change Sample name SERUM to Serum
#   + Create SD column and remove MEAN+SD column
#   + For each sample type
#     - Round time to closest value out of nominal times
#   + For each sample type and nominal time
#     - Find the geometric mean of the MEAN and SD values (Geometric mean is 
#       used as data was digitised on a log-scale, i.e. clicking a pixel above
#       the point leads to greater error than a pixel below)
# * Combine datasets by row
  figdf <- map_dfr(c("fig1", "fig2"), ~ .x %>%
    paste("olsendig", sep = "_") %>%
    get() %>%
    rename(Time = x) %>%
    mutate_at(vars(matches("MEAN")), ~ .x * 1000) %>%
    pivot_longer(cols = matches("MEAN"), names_to = "Sample_Metric") %>%
    separate(Sample_Metric, c("Sample", "Metric"), sep = "_") %>%
    pivot_wider(id_cols = c("Sample", "Time"), names_from = "Metric") %>%
    mutate(Sample = if_else(Sample == "SERUM", "Serum", Sample), 
      SD = abs(MEANPLUSSD - MEAN)) %>%
    select(-MEANPLUSSD) %>%
    group_by(Sample) %>%
    mutate(Time = map_dbl(Time, ~ .x %>%
      magrittr::subtract(unique(tabdf$Time)) %>%
      abs() %>%
      {which(. == min(.))} %>%
      {magrittr::extract2(unique(tabdf$Time), .)}
    )) %>%
    group_by(Sample, Time) %>%
    summarise_all(~ exp(mean(log(.x))))
  )
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Determine if table data is the data that was plotted in the figures
# Combine the two datasets and calculate the relative difference in values
# * Combine the two dataset with an identifier column to differentiate source
# * Reshape MEAN and SD columns two columns: column name (Metric) and value
# * Combine Data source and Metric into one column
# * Transform each unique combindation of Data source and Metric into its own column
# * Filter out times for each sample that were zero for the table data
# * Calculate the relative difference between the two datasets
  combdf <- tabdf %>%
    add_column(Data = "tab") %>%
    bind_rows(add_column(figdf, Data = "fig")) %>%
    pivot_longer(cols = c("MEAN", "SD"), names_to = "Metric") %>%
    unite(col = "Data_Metric", Data, Metric, sep = "_") %>%
    pivot_wider(id_cols = c("Sample", "Time"), names_from = "Data_Metric") %>%
    filter(tab_MEAN != 0) %>%
    mutate(dMEAN = tab_MEAN - fig_MEAN, dSD = tab_SD - fig_SD) %>%
    mutate(dMEANrel = 100*dMEAN/tab_MEAN, dSDrel = 100*dSD/tab_SD) %>%
    select(Sample, Time, dMEANrel, dSDrel)
# Summarise the relative difference for each Sample and each Time (mean)
  combdf %>%
    select(-Time) %>%
    group_by(Sample) %>%
    summarise_all(mean, na.rm = TRUE)
  combdf %>%
    select(-Sample) %>%
    group_by(Time) %>%
    summarise_all(mean, na.rm = TRUE)
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# The table data is what was provided in the plots
# This dataset is the cleanest and has the least error so will be used for
# comparison to model predictions
  write_csv(tabdf, "OlsenTable1.csv")
  