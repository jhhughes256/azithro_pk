# Azithromycin Population Pharmacokinetic Model
# ------------------------------------------------------------------------------
# Model sourced from:
#   Songmao Zheng, Peter Matzneller, Markus Zeitlinger, Stephan Schmidt, 
#   Development of a Population Pharmacokinetic Model Characterizing the Tissue
#   Distribution of Azithromycin in Healthy Subjects, Antimicrobial Agents and 
#   Chemotherapy Oct 2014, 58 (11) 6675-6684; DOI: 10.1128/AAC.02904-14
# This script aims to replicate Figure 3 from that manuscript.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load packages
  library(tidyverse)
  library(mrgsolve)
  library(scales)

# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.margin = unit(c(1, 1.5, 1, 1), "lines"))
  
# Define colourblind palette index
  cbPal <- data.frame(
    blue = "#0093D0",  # pfizerblue
    lblue = "#A6CEE3",
    green = "#33A02C",
    lgreen = "#B2DF8A",
    pink = "#FB9A99",
    red = "#E31A1C",
    orange = "#FF7F00",
    lorange = "#FDBF6F",
    purple = "#6A3D9A",
    lpurple = "#CAB2D6",
    yellow = "#FFFF99",
    brown = "#B15928",
    stringsAsFactors = F
  )
  
# Source model code and compile (compiled model object: `mod` )
  source("evaluation/01_InitialTest/model.R")
  
# Create population of individuals
  idata <- omat(mod, make = TRUE) %>%
    MASS::mvrnorm(n = 1000, mu = rep(0, nrow(.)), Sigma = .) %>%
    as_tibble() %>%
    set_names(paste0("ETA", 1:ncol(.))) %>%
    mutate(ID = 1:1000)

# Run a trial simulation
  simdf <- mod %>% 
    idata_set(idata) %>%
    ev(amt = 500, ii = 24, addl = 2) %>% 
    mrgsim(end = 220, delta = 0.5) %>%
    as_tibble() %>%
    select(time, CPLAST, CMUSC, CSUBC, CPMLT) %>%
    pivot_longer(cols = contains("C"), names_to = "DVID", values_to = "DV") %>%
    group_by(time, DVID) %>%
    summarise(DVmean = mean(DV), DVlo = quantile(), DVhi) %>%
    ungroup() %>%
    mutate(DVIDf = factor(DVID, labels = c("Muscle (unb.)", "Plasma (tot.)", 
      "PML (tot.)", "Subcutis (unb.)")))
  
# Plot results of trial simulation
  p <- NULL
  p <- ggplot(simdf)
  p <- p + geom_line(aes(x = time, y = DVmean, colour = DVIDf), size = 1)
  p <- p + geom_hline(yintercept = 2000, linetype = "dashed", 
    colour = cbPal$brown, size = 1)
  p <- p + geom_text(x = 204, y = log10(2000)-0.3, label = "~MIC[90]",
    parse = TRUE)
  p <- p + scale_x_continuous("Time (hour)", breaks = 0:9*24)
  p <- p + scale_y_log10("Concentration (ng/mL)", 
    breaks = c(1,  10, 100, 1000, 10000, 100000), labels = comma)
  p <- p + coord_cartesian(xlim = c(0, 216), ylim = c(1, 100000), expand = FALSE)
  p <- p + scale_colour_manual("Tissue", 
    values = with(cbPal, c(green, orange, blue, red)))
  p <- p + theme(legend.position = "bottom")
  p  
  
