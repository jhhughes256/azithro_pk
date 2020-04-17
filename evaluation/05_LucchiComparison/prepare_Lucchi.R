# Prepare digitised azithromycin data from Lucchi paper
# -----------------------------------------------------------------------------
# This script is designed to prepare digitised plot data to compare with 
# predictions from the Zheng S et al. mrgsolve model.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data digitised from:
#   Lucchi, M., et al. (2008). "Pharmacokinetics of azithromycin in serum, 
#   bronchial washings, alveolar macrophages and lung tissue following a single 
#   oral dose of extended or immediate release formulations of azithromycin." 
#   J Antimicrob Chemother 61(4): 884-891.
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
  fig1_lucchidig <- read_csv("Digitised_LucchiFigure1.csv", 
    col_types = cols(.default = col_double()))
  fig2_lucchidig <- read_csv("Digitised_LucchiFigure2.csv", 
    col_types = cols(.default = col_double()))
  fig3_lucchidig <- read_csv("Digitised_LucchiFigure3.csv", 
    col_types = cols(.default = col_double()))
  fig4_lucchidig <- read_csv("Digitised_LucchiFigure4.csv",
    col_types = cols(.default = col_double()))

# Define nominal times
  nom_times <- c(0, 2, 4, 8, 12, 16, 24, 48, 72)
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare digitised data to be compared
# Digitised data (units: mg/L)
# * For each dataset
#   + Rename x to Time and convert values to ng/mL
#   + Create unique row ID to prevent reshape errors
#   + Reshape all columns except Time into two columns: column name and value
#   + Split column names into Form (IR, ER, etc.) and Metric (MEAN etc.)
#   + Transform each unique Metric into its own column
#   + Create Sample variable according to source data
#   + Create SD column and remove MEAN+SD and UID column
#   + For each sample type
#     - Round time to closest value out of nominal times
#   + Remove pre-dose concentrations (nominal time == 0)
#   + For each sample type and nominal time
#     - Find the mean of the MEAN and SD values
# * Combine datasets by row
  figdf <- map_dfr(c("fig1", "fig2", "fig3", "fig4"), ~ .x %>%
    paste("lucchidig", sep = "_") %>%
    get() %>%
    rename(Time = x) %>%
    mutate_at(vars(matches("MEAN")), ~ .x * 1000) %>%
    mutate(UID = 1:nrow(.)) %>%
    pivot_longer(cols = matches("MEAN"), names_to = "Form_Metric", values_to = "value") %>%
    separate(Form_Metric, c("Form", "Metric"), sep = "_") %>%
    pivot_wider(id_cols = c("Time", "Form", "UID"), names_from = "Metric", values_from = "value") %>%
    add_column(Sample = case_when(
      .x == "fig1" ~ "Serum",
      .x == "fig2" ~ "ELF",
      .x == "fig3" ~ "AM",
      .x == "fig4" ~ "Lung"
    ), .before = "Form") %>%
    mutate(SD = abs(MEANPLUSSD - MEAN)) %>%
    select(-MEANPLUSSD, -UID) %>%
    group_by(Sample, Form) %>%
    mutate(Time = map_dbl(Time, ~ .x %>%
      magrittr::subtract(nom_times) %>%
      abs() %>%
      {which(. == min(.))} %>%
      {magrittr::extract2(unique(nom_times), .)}
    )) %>%
    filter(Time != 0) %>%
    group_by(Sample, Form, Time) %>%
    summarise_all(mean)
  )

# Save digitised data
  write_csv(figdf, "LucchiFigure123.csv")
  