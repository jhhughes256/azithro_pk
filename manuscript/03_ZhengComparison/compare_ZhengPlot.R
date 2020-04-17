# Azithromycin Population Pharmacokinetic Model
# ------------------------------------------------------------------------------
# Model sourced from:
#   Songmao Zheng, Peter Matzneller, Markus Zeitlinger, Stephan Schmidt, 
#   Development of a Population Pharmacokinetic Model Characterizing the Tissue
#   Distribution of Azithromycin in Healthy Subjects, Antimicrobial Agents and 
#   Chemotherapy Oct 2014, 58 (11) 6675-6684; DOI: 10.1128/AAC.02904-14
# This script aims to replicate Figure 3 from that manuscript.
# Data digitised using open source digitising software Engauge v12.1
# Available at: https://github.com/markummitchell/engauge-digitizer/releases
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Regimen from Zheng et al. is: 500mg D1-D3
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare work environment
# Clear objects from R environment
  rm(list = ls(all = TRUE))

# Load packages
  library(tidyverse)
  library(mrgsolve)
  library(scales)

# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.margin = unit(c(1, 1.5, 1, 1), "lines"))
  
# Source external files
# * Compile mrgsolve model
  source("model.R")  # `mrgsolve` model object: mod 

# * Read in digitised data
  line_zhengdig <- read_csv("DigitisedLines_Zhengplot.csv", 
    col_types = cols(.default = col_double()))
  plas_zhengdig <- read_csv("DigitisedObs_PLAS_Zhengplot.csv", 
    col_types = cols(.default = col_double()))
  musc_zhengdig <- read_csv("DigitisedObs_MUSC_Zhengplot.csv", 
    col_types = cols(.default = col_double()))
  subc_zhengdig <- read_csv("DigitisedObs_SUBC_Zhengplot.csv", 
    col_types = cols(.default = col_double()))
  pmlc_zhengdig <- read_csv("DigitisedObs_PML_Zhengplot.csv", 
    col_types = cols(.default = col_double()))
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare digitised data for plotting
# Zheng et al. line data
  linedf <- line_zhengdig %>%
    rename(time = x) %>%
    pivot_longer(cols = contains("LINE"), names_to = "DVID", values_to = "DV") %>%
    mutate(DVIDf = factor(DVID, labels = c("Muscle (unb.)", "Plasma (tot.)", 
      "PML (tot.)", "Subcutis (unb.)")))
  
# Zheng et al. error bar data
  obsdf <- map_dfr(c("PLAS", "MUSC", "SUBC", "PMLC"), function(cmt) {
    get(paste0(tolower(cmt), "_zhengdig")) %>%
      rename(OBS = paste0(cmt, "_OBS"), SE = paste0(cmt, "_SE")) %>%
      mutate(DVID = cmt)}) %>% 
    rename(time = x) %>%
    mutate(DVIDf = factor(DVID, labels = c("Muscle (unb.)", "Plasma (tot.)", 
      "PML (tot.)", "Subcutis (unb.)")))
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Create Zheng-style Plot
# Uses total plasma and PML concentrations instead of unbound/unionised
# Run a trial simulation
  simdf <- mod %>% 
    ev(amt = 500, ii = 24, addl = 2) %>% 
    mrgsim(end = 220, delta = 0.5) %>%
    as_tibble() %>%
    select(time, CPLAST, CMUSC, CSUBC, CPMLT) %>%
    pivot_longer(cols = contains("C"), names_to = "DVID", values_to = "DV") %>%
    mutate(DVIDf = factor(DVID, labels = c("Muscle (unb.)", "Plasma (tot.)", 
      "PML (tot.)", "Subcutis (unb.)")))
  
# Plot results of trial simulation
  p <- NULL
  p <- ggplot(simdf)
  p <- p + geom_line(aes(x = time, y = DV, colour = DVIDf, linetype = "Model"), 
    size = 1)
  p <- p + geom_line(aes(x = time, y = DV, colour = DVIDf, linetype = "Zheng"), 
    size = 1, data = linedf)
  p <- p + geom_point(aes(x = time, y = OBS, colour = DVIDf), 
    size = 2, shape = 1, data = obsdf)
  p <- p + scale_x_continuous("Time (hour)", breaks = 0:9*24)
  p <- p + scale_y_log10("Concentration (ng/mL)", 
    breaks = c(1,  10, 100, 1000, 10000, 100000), labels = comma)
  p <- p + coord_cartesian(xlim = c(0, 216), ylim = c(1, 100000), expand = FALSE)
  p <- p + scale_colour_manual("Tissue", 
    values = c("darkgreen", "orange", "blue", "red"))
  p <- p + scale_linetype_manual("", 
    values = c("solid", "dotted"))
  p <- p + theme(legend.position = "bottom")
  p  
  
