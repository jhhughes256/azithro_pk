# Azithromycin Population Pharmacokinetic Model
# ------------------------------------------------------------------------------
# Regimen from Baldwin et al. 500 mg single dose
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Data digitised from:
#   Baldwin, D. R., et al. (1990). "Azithromycin concentrations at the sites of 
#   pulmonary infection." Eur Respir J 3(8): 886-890.
# Data digitised using open source digitising software Engauge v12.1
# Available at: https://github.com/markummitchell/engauge-digitizer/releases
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
  source("../../model.R")  # `mrgsolve` model object: mod 

# * Read in prepared digitised data
  baldwindf <- read_csv("BaldwinFigure2.csv", 
    col_types = cols(.default = col_double(), Sample = col_character())) %>%
    mutate(DVIDf = factor(Sample, labels = c("AMs", "ELF", "Plasma")))
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Create Zheng-style Plot
# Uses total plasma and PML concentrations instead of unbound/unionised
# Run a trial simulation
  simdf <- mod %>% 
    ev(amt = 500) %>% 
    mrgsim(end = 96, delta = 0.5) %>%
    as_tibble() %>%
    select(time, CPLAST, CMUSC, CSUBC, CPMLT, CPMLLT) %>%
    pivot_longer(cols = contains("C"), names_to = "DVID", values_to = "DV") %>%
    mutate(DVIDf = factor(DVID, labels = c("Muscle", "Plasma", "PML Lysosome", 
      "PML", "Subcutis")))
  
# Plot results of trial simulation
  p <- NULL
  p <- ggplot(simdf)
  p <- p + geom_line(aes(x = time, y = DV, colour = DVIDf), size = 1)
  p <- p + geom_point(aes(x = Time, y = MEAN, colour = DVIDf), 
    size = 2, shape = 1, data = baldwindf)
  p <- p + geom_errorbar(aes(x = Time, ymin = MEAN-SD, ymax = MEAN+SD, 
    colour = DVIDf), size = 0.5, width = 1, data = baldwindf)
  p <- p + scale_x_continuous("Time (hour)", breaks = 0:8*12)
  p <- p + scale_y_log10("Concentration (ng/mL)", 
    breaks = c(1,  10, 100, 1000, 10000, 100000, 1000000), labels = comma)
  p <- p + coord_cartesian(xlim = c(0, 96), ylim = c(1, 1000000))
  p <- p + scale_colour_manual("Tissue",
    values = c("purple", "goldenrod", "blue", "brown", "red", "pink", "black"))
  p <- p + scale_linetype_manual("", 
    values = c("solid", "dotted"))
  p <- p + theme(legend.position = "bottom")
  p  
  
