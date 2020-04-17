# Prepare data from digitised sources
# -----------------------------------------------------------------------------
# Digitised data sourced from the following publications:
#   Ballow, CH, et al. (1998). "Pharmacokinetics of Oral Azithromycin in...
#   Danesi, R, et al. (2003). "Comparative distribution of azithromycin in...
#   Di Paolo, A, et al. (2002). "PHARMACOKINETICS OF AZITHROMYCIN IN LUNG...
#   Liu, P et al. (2007). "Comparative pharmacokinetics of azithromycin in...
#   Lucchi, M, et al. (2008). "Pharmacokinetics of azithromycin in serum...
#   Olsen, KM, et al. (1996). "Intrapulmonary pharmacokinetics of...
#   Sampson, MR, et al. (2014). "Population Pharmacokinetics of Azithromycin...
#   Zheng, S, et al. (2014). "Development of a population pharmacokinetic...
# Figure data digitised using open source digitising software Engauge v12.1
# Available at: https://github.com/markummitchell/engauge-digitizer/releases
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Prepare work environment
# Clear objects from R environment
  rm(list = ls(all = TRUE))

# Load packages
  library(tidyverse)
  
# Source external files
# * Read in Ballow digitised data and add metadata
  balldf <- read_csv("BallowFigure1.csv", 
      col_types = cols(.default = col_double(), Sample = col_character())) %>%
    group_by(Time, Sample) %>%
    summarise(Mean = mean(Obs), SD = sd(Obs)) %>%
    rename(DV = Mean) %>%
    mutate(Source = "Ballow et al.", DVTYPE = "Obs", SDTYPE = "SD") %>%
    mutate(Form = "IR", Regimen = "D1_500mg_D2-D5_250mg") %>%
    mutate(Dose = if_else(Time < 24, 500, 250))
  
# * Read in Danesi digitised data and add metadata
  danedf <- read_csv("DanesiFigure1.csv", 
      col_types = cols(.default = col_double(), Sample = col_character())) %>%
    rename(DV = Mean) %>%
    mutate(Source = "Danesi et al.", DVTYPE = "Mean", SDTYPE = "SD") %>%
    mutate(Form = "IR", Regimen = paste0("D1-D3_", Dose, "mg"))
  
# * Read in Di Paolo digitised data and add metadata
  dipadf <- read_csv("DiPaoloFigure1.csv", 
      col_types = cols(.default = col_double(), Sample = col_character())) %>%
    rename(DV = Mean) %>%
    mutate(Source = "Di Paolo et al.", DVTYPE = "Mean", SDTYPE = "SD") %>%
    mutate(Form = "IR", Regimen = paste0("D1-D3_", Dose, "mg"))
  
# * Read in Liu digitised data and add metadata
  liupdf_path <- "A0661112FigureS123.csv"
  if ("A0661112_summary.csv" %in% dir()) { liupdf_path <- "A0661112_summary.csv" }
  liupdf <- read_csv(liupdf_path, 
      col_types = cols(.default = col_double(), 
        Sample = col_character(), Form = col_character())) %>%
    rename(DV = Mean) %>%
    mutate(Source = "Liu et al.", DVTYPE = "Mean", SDTYPE = "SD") %>%
    mutate(Dose = if_else(Form == "SR", 2000, 500)) %>%
    mutate(Regimen = paste0(if_else(Form == "SR", "D1_", "D1-D3_"), Dose, "mg"))
  
# * Read in Lucchi digitised data and add metadata
  luccdf <- read_csv("LucchiFigure123.csv", 
      col_types = cols(.default = col_double(), 
        Sample = col_character(), Form = col_character())) %>%
    rename(DV = MEAN) %>%
    mutate(Source = "Lucchi et al.", DVTYPE = "Mean", SDTYPE = "SD") %>%
    mutate(Form = if_else(Form == "ER", "SR", Form)) %>%
    mutate(Dose = if_else(Form == "SR", 2000, 500)) %>%
    mutate(Regimen = paste0("D1_", Dose, "mg"))
  
# * Read in Olsen digitised data and add metadata
  olsedf <- read_csv("OlsenTable1.csv", 
      col_types = cols(.default = col_double(), Sample = col_character())) %>%
    rename(DV = MEAN) %>%
    mutate(Source = "Olsen et al.", DVTYPE = "Mean", SDTYPE = "SD") %>%
    mutate(Form = "IR", Regimen = "D1_500mg_D2-D5_250mg") %>%
    mutate(Dose = if_else(Time < 24, 500, 250))
  
# * Read in Sampson digitised data and add metadata
  sampdf <- read_csv("SampsonFigure1AB.csv", 
      col_types = cols(.default = col_double(), Sample = col_character())) %>%
    rename(DV = Mean, SD = SE) %>%
    mutate(Source = "Sampson et al.", DVTYPE = "Mean", SDTYPE = "SE") %>%
    mutate(Form = "IR", Regimen = paste0("D1_", Dose, "mg"))
  
# * Read in Zheng digitised data and add metadata
# Renaming to Matzneller to align with manuscript
  zhendf <- read_csv("ZhengFigure3.csv", 
      col_types = cols(.default = col_double(), Sample = col_character())) %>%
    rename(DV = MEAN, SD = SE) %>%
    mutate(Source = "Matzneller et al.", DVTYPE = "Mean", SDTYPE = "SE") %>%
    mutate(Form = "IR", Dose = 500, Regimen = "D1-D3_500mg")
  
# Bind data together and standardise/filter Sample column
  outdf <- bind_rows(balldf, danedf, dipadf, liupdf, luccdf, olsedf, sampdf, zhendf) %>%
    select(Source, Form, Regimen, Dose, Sample, Time, DV, SD, DVTYPE, SDTYPE) %>%
    mutate(Sample = case_when(
      Sample == "MUSC" ~ "Muscle",
      Sample %in% c("PLAS", "Plasma", "Serum") ~ "Plasma (tot.)",
      Sample %in% c("PMN", "PML", "PMLC") ~ "PML",
      Sample == "SUBC" ~ "Subcutis",
      TRUE ~ Sample)) %>%
    filter(Sample %in% c("AM", "Lung", "MNL", "Muscle", "PBM", "Plasma (tot.)", "PML", "Subcutis"))
    
# Write data to file
  write_csv(outdf, "AZDigitised_PK_10APR2020.csv")
  