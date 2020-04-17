# Azithromycin Population Pharmacokinetic Model
# ------------------------------------------------------------------------------
# This script aims to validate mrgsolve model code against the NONMEM code 
# sourced from:
#   Songmao Zheng, Peter Matzneller, Markus Zeitlinger, Stephan Schmidt, 
#   Development of a Population Pharmacokinetic Model Characterizing the Tissue
#   Distribution of Azithromycin in Healthy Subjects, Antimicrobial Agents and 
#   Chemotherapy Oct 2014, 58 (11) 6675-6684; DOI: 10.1128/AAC.02904-14
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# The regimen used to simulate this data is from:
#   Gautret et al. Hydroxychloroquine and azithromycin as a treatment of 
#   COVID-19: results of an open-label non-randomized clinical trial, 
#   International Journal of Antimicrobial Agents, Mar 2020; 
#   DOI: 10.1016/j/ijantimicag.2020.105949
# The regimen for azithromycin is: 500mg D1, 250mg D2-5
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare work environment
# Remove all current objects in the workspace
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
  source("evaluation/02_NONMEMValidation/model.R")  # `mrgsolve` model object: mod 

# * Read in NONMEM input
  obsdf <- read_csv("evaluation/02_NONMEMValidation/AZPop_30MAR2020.csv", na = ".", 
    col_types = cols(.default = col_double(), C = col_character()))
  nid <- length(unique(obsdf$ID))
  
# * Read in NONMEM simulation output
  simdf <- list.files("evaluation/02_NONMEMValidation") %>% # List files in working directory
    str_which(".fit") %>% # Identify position of file with ".fit" extension
    {list.files()[.]} %>% # Identify filename of file in x position
    read_table(na = c("NA","***********", "1.#INFE+00"), skip = 1,
      col_types = cols(.default = col_double())) %>%
    filter(EVID == 0) %>%  # subset for PK data
    mutate(PROJ = "B793") %>%  # redefine PROJ number
    select(PROJ, everything(), -MDV, -EVID)  # remove MDV EVID
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulation and Validation of PRED values
# Create input dataset
# * Remove C column
# * Filter CMT values (DV reported via individual columns)
  mrgindf1 <- obsdf %>%
    select(-C) %>%
    filter(CMT %in% c(1, 2))
  
# Simulate PRED using `mrgsolve` model
# * Identify input dataset and variables to carry through to the output
# * Simulate and filter out dose observations
  mrgpreddf <- mod %>%
    data_set(mrgindf1) %>%
    carry_out(amt, evid, rate, cmt) %>%
    mrgsim_df() %>%
    filter(evid == 0)
  
# Validation of `mrgsolve` model output
# Numeric difference and relative difference
# * For each combination of CMT and mrgsolve PRED variable name
# * Subset NONMEM output for CMT
# * Convert PRED from that subset to ng/mL and add to mrgsolve output
# * Determine range of actual and relative differences
  map2(c(2, 3, 4, 9), c("CPLAS", "CMUSC", "CSUBC", "CPMLC"), ~ simdf %>%
    filter(CMT == .x) %>%
    {mutate(mrgpreddf, NMPRED = 1000*.$PRED)} %>%
    select(ID, TIME, !!as.name(.y), NMPRED) %>%
    with(list(
      numeric = range(NMPRED - get(.y)) %>% signif(3),
      relative = range(100*(NMPRED - get(.y))/(NMPRED+1e-26)) %>%
        signif(3) %>%
        paste("%")
    ))
  ) %>% set_names(c("CPLAS", "CMUSC", "CSUBC", "CPMLC"))
  
# Visual evaluation
# Compare PRED over the simulated time-course
# * Define labels for compartment concentrations
# * Set up pdf graphics device
# * Pivot mrgsolve dataset so that PRED is all in one column
# * Create factor column for compartment for labelling
# * Plot PRED for each of the four compartments for both models
# * Turn off graphics device
  DV_labs <- c("Plasma", "Muscle", "Subcutis", "PML cytosol")
  pdf("evaluation/02_NONMEMValidation/PRED_vs_TIME.pdf", width = 8, height = 10)
  p <- mrgpreddf %>%
    filter(ID == 1) %>%
    pivot_longer(cols = c("CPLAS", "CMUSC", "CSUBC", "CPMLC"), 
      names_to = "DVID", values_to = "PRED") %>%
    mutate(DVIDf = factor(DVID, labels = DV_labs)) %>%
    ggplot()
  p <- p + geom_line(aes(x = TIME, y = PRED, colour = DVIDf, linetype = "mrgsolve"), 
    size = 1, alpha = 0.5)
  p <- p + geom_line(aes(x = TIME, y = 1000*PRED, colour = DVIDf, linetype = "NONMEM"), 
    size = 1, alpha = 0.5, data = simdf %>%
      filter(ID == 1) %>%
      mutate(DVIDf = factor(CMT, labels = DV_labs))
  )
  p <- p + scale_x_continuous("Time (hour)", breaks = 0:5*24)
  p <- p + scale_y_log10("Concentration (ng/mL)", 
    breaks = c(0.1, 0.3, 1, 3, 10, 30, 100, 300))
  # p <- p + coord_cartesian(xlim = c(0, 120))
  p <- p + scale_colour_manual("Tissue", 
    values = c("darkgreen", "goldenrod", "blue", "red"))
  p <- p + scale_linetype_manual("Model", 
    values = c("solid", "dotted"))
  p <- p + theme(legend.position = "bottom")
  p  
  dev.off()
   
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Simulation and Validation of IPRED values
# Create input dataset
# * Extract ETA values simulated by NONMEM
# * Rename to corresponding ETAs in mrgsolve model
# * Merge input dataset with ETA dataset
# * Remove C column
# * Filter CMT values (DV reported via individual columns)
  mrgindf2 <- simdf %>%
    distinct(ID, ETA2, ETA3, ETA4, ETA6, ETA7, ETA8, ET15) %>%
    set_names(c("ID", paste0("ETA", c(2:4, 1, 5:7)))) %>%
    {full_join(obsdf, ., by = "ID")} %>%
    select(-C) %>%
    filter(CMT %in% c(1, 2))
  
# Simulate IPRED using `mrgsolve` model
# * Identify input dataset and variables to carry through to the output
# * Simulate and filter out dose observations
  mrgipredf <- mod %>%
    data_set(mrgindf2) %>%
    carry_out(amt, evid, rate, cmt) %>%
    mrgsim_df() %>%
    filter(evid == 0)
  
# Validation of `mrgsolve` model output
# Numeric difference and relative difference
# * For each combination of CMT and mrgsolve IPRED variable name
# * Subset NONMEM output for CMT
# * Convert IPRED from that subset to ng/mL and add to mrgsolve output
# * Determine range of actual and relative differences
  map2(c(2, 3, 4, 9), c("CPLAS", "CMUSC", "CSUBC", "CPMLC"), ~ simdf %>%
    filter(CMT == .x) %>%
    {mutate(mrgipredf, NMIPRED = 1000*.$IPRED)} %>%
    select(ID, TIME, !!as.name(.y), NMIPRED) %>%
    with(list(
      numeric = range(NMIPRED - get(.y)) %>% signif(3),
      relative = range(100*(NMIPRED - get(.y))/(NMIPRED+1e-26)) %>%
        signif(3) %>%
        paste("%")
    ))
  ) %>% set_names(c("CPLAS", "CMUSC", "CSUBC", "CPMLC"))
  
# Visual evaluation
# Compare mean IPRED over the simulated time-course
# * Set up pdf graphics device
# * Pivot mrgsolve dataset so that PRED is all in one column
# * Create factor column for compartment for labelling
# * Plot mean IPRED for each of the four compartments for both models
# * Turn off graphics device
  pdf("evaluation/02_NONMEMValidation/meanIPRED_vs_TIME.pdf", width = 8, height = 10)
  p <- mrgipredf %>%
    pivot_longer(cols = c("CPLAS", "CMUSC", "CSUBC", "CPMLC"), 
      names_to = "DVID", values_to = "IPRED") %>%
    mutate(DVIDf = factor(DVID, labels = DV_labs)) %>%
    ggplot()
  p <- p + stat_summary(aes(x = TIME, y = IPRED, colour = DVIDf, 
    linetype = "mrgsolve"), size = 1, alpha = 0.5, geom = "line", 
    fun.y = "mean", fun.args = list(na.rm = TRUE))
  p <- p + stat_summary(aes(x = TIME, y = 1000*IPRED, colour = DVIDf, 
    linetype = "NONMEM"), size = 1, alpha = 0.5, geom = "line", 
    fun.y = "mean", fun.args = list(na.rm = TRUE),
    data = simdf %>%
      mutate(DVIDf = factor(CMT, labels = DV_labs))
  )
  p <- p + scale_x_continuous("Time (hour)", breaks = 0:5*24)
  p <- p + scale_y_log10("Concentration (ng/mL)", 
    breaks = c(0.1, 0.3, 1, 3, 10, 30, 100, 300))
  # p <- p + coord_cartesian(xlim = c(0, 120))
  p <- p + scale_colour_manual("Tissue", 
    values = c("darkgreen", "goldenrod", "blue", "red"))
  p <- p + scale_linetype_manual("Model", 
    values = c("solid", "dotted"))
  p <- p + theme(legend.position = "bottom")
  p
  dev.off()
  
# Visual evaluation
# Compare median IPRED and confidence intervals for each compartment
# * Set up pdf graphics device
# * For each compartment and IPRED name
# * Subset NONMEM output for CMT
# * Convert IPRED from that subset to ng/mL and add to mrgsolve output
# * Pivot mrgsolve dataset so that IPRED is all in one column
# * Calculate median and 90% confidence intervals at each TIME for each MODEL
# * Create factor column for model for labelling
# * Plot median IPRED and 90% CI over time
# * Turn off graphics device
  pdf("evaluation/02_NONMEMValidation/ciIPRED_vs_TIME.pdf", width = 8, height = 10)
  map2(c(2, 3, 4, 9), c("CPLAS", "CMUSC", "CSUBC", "CPMLC"), function(cmt, par) {
    p <- simdf %>%
      filter(CMT == cmt) %>%
      {mutate(mrgipredf, NMIPRED = 1000*.$IPRED)} %>%
      select(ID, TIME, !!as.name(par), NMIPRED) %>%
      pivot_longer(cols = c(!!as.name(par), NMIPRED), 
        names_to = "MODEL", values_to = "IPRED") %>%
      group_by(TIME, MODEL) %>%
      summarise(median = median(IPRED), 
        ci90lo = quantile(IPRED, prob = 0.05), 
        ci90hi = quantile(IPRED, prob = 0.95)) %>%
      mutate(MODELf = factor(MODEL, labels = c("mrgsolve", "NONMEM"))) %>%
      ggplot()
    p <- p + geom_line(aes(x = TIME, y = median, colour = MODELf, linetype = "Median"), 
      size = 1, alpha = 0.5)
    p <- p + geom_line(aes(x = TIME, y = ci90lo, colour = MODELf, linetype = "90% CI"), 
      size = 1, alpha = 0.5)
    p <- p + geom_line(aes(x = TIME, y = ci90hi, colour = MODELf, linetype = "90% CI"), 
      size = 1, alpha = 0.5)
    p <- p + scale_x_continuous("Time (hour)", breaks = 0:5*24)
    p <- p + scale_y_log10(paste(par, "Concentration (ng/mL)"))
    # p <- p + coord_cartesian(xlim = c(0, 120))
    p <- p + scale_colour_manual("Model", 
      values = c("red", "blue"))
    p <- p + scale_linetype_manual("", values = c("dashed", "solid"))
    p <- p + guides(colour = guide_legend(override.aes = list(alpha = 1))) 
    p <- p + guides(linetype = guide_legend(override.aes = list(alpha = 1))) 
    p <- p + theme(legend.position = "bottom")
    p
  })
  dev.off()
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Simulation and Evaluation of OMEGA distribution
# Create input dataset
# * Using mrgsolve OMEGA values sample ETA from normal distribution
# * Rename to corresponding ETAs in mrgsolve model
# * Add ID column and merge input dataset with ETA dataset
# * Remove C column
# * Filter CMT values (DV reported via individual columns)
  mrgindf3 <- mod %>%  
    omat(make = TRUE) %>%
    {MASS::mvrnorm(n = nid, mu = rep(0, times = nrow(.)), Sigma = .)} %>%
    as_tibble() %>%
    set_names(paste0("ETA", 1:7)) %>%
    add_column(ID = 1:nid, .before = "ETA1") %>%
    {full_join(obsdf, ., by = "ID")} %>%
    select(-C) %>%
    filter(CMT %in% c(1, 2))
    
# Simulate individual parameters using `mrgsolve` model
# * Identify input dataset and variables to carry through to the output
# * Simulate and filter out dose observations
  mrgomegdf <- mod %>%
    data_set(mrgindf3) %>%
    carry_out(amt, evid, rate, cmt) %>%
    mrgsim_df() %>%
    filter(evid == 0)
  
# Evaluation of `mrgsolve` model output
# Visual comparison of parameter distributions
# * Define corresponding parameter names for NONMEM and mrgsolve output
# * Set up pdf graphics device
# * For each combination of parameter names
# * Plot histograms of NONMEM and mrgsolve parameter values
# * Turn off graphics device
  nmpar <- c("TVLAG", "CL", "V2", "KIN", "KPMU", "KPSU", "KPLE")
  mrgpar <- c("TLAG", "CLF", "VCF", "KIN", "DFMUSC", "DFSUBC", "DFPMLC")
  pdf("evaluation/02_NONMEMValidation/THETAi_hist.pdf", width = 4, height = 5)
  map2(nmpar, mrgpar, function(nm, mrg) {
    p <- ggplot()
    p <- p + geom_histogram(aes_string(x = nm), fill = "red", alpha = 0.5,  
      data = distinct(simdf, ID, .keep_all = TRUE))
    p <- p + geom_histogram(aes_string(x = mrg), fill = "blue", alpha = 0.5,
      data = distinct(mrgomegdf, ID, .keep_all = TRUE))
    p <- p + scale_x_continuous(paste(mrg, "n =", nid))
    p <- p + scale_y_continuous("Count")
    p
  })
  dev.off()
  pdf("evaluation/02_NONMEMValidation/THETAi_dens.pdf", width = 4, height = 5)
  map2(nmpar, mrgpar, function(nm, mrg) {
    p <- ggplot()
    p <- p + geom_histogram(aes_string(x = nm), fill = "red", alpha = 0.5,  
      stat = "density", data = distinct(simdf, ID, .keep_all = TRUE))
    p <- p + geom_histogram(aes_string(x = mrg), fill = "blue", alpha = 0.5,
      stat = "density", data = distinct(mrgomegdf, ID, .keep_all = TRUE))
    p <- p + scale_x_continuous(paste(mrg, "n =", nid))
    p <- p + scale_y_continuous("Density")
    p
  })
  dev.off()
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Simulation and Evaluation of SIGMA distribution (via residuals)
# Create residual column for mrgsolve output
  mrgresdf <- mutate(mrgipredf, 
    IRESPLAS = CPLAS - DVPLAS, IRESMUSC = CMUSC - DVMUSC,
    IRESSUBC = CSUBC - DVSUBC, IRESPMLC = CPMLC - DVPMLC
  )
  
# Evaluation of `mrgsolve` model output
# Visual comparison of residual distributions
# * Set up pdf graphics device
# * For each combination of compartment number and name
# * Plot histograms of NONMEM and mrgsolve residuals
# * Turn off graphics device
  pdf("evaluation/02_NONMEMValidation/IRES_hist.pdf", width = 4, height = 5)
  map2(c(2, 3, 4, 9), c("PLAS", "MUSC", "SUBC", "PMLC"), function(cmt, par) {
    p <- ggplot()
    p <- p + geom_histogram(aes(x = 1000*IRES), fill = "red", alpha = 0.5,  
      data = filter(simdf, CMT == cmt))
    p <- p + geom_histogram(aes_string(x = paste0("IRES", par)), fill = "blue", alpha = 0.5,
      data = mrgresdf)
    p <- p + scale_x_continuous(paste(par, "Residuals"))
    p <- p + scale_y_continuous("Count")
    p
  })
  pdf("evaluation/02_NONMEMValidation/IRES_dens.pdf", width = 4, height = 5)
  map2(c(2, 3, 4, 9), c("PLAS", "MUSC", "SUBC", "PMLC"), function(cmt, par) {
    p <- ggplot()
    p <- p + geom_histogram(aes(x = 1000*IRES), fill = "red", alpha = 0.5,  
      stat = "density", data = filter(simdf, CMT == cmt))
    p <- p + geom_histogram(aes_string(x = paste0("IRES", par)), fill = "blue", alpha = 0.5,
      stat = "density", data = mrgresdf)
    p <- p + scale_x_continuous(paste(par, "Residuals"))
    p <- p + scale_y_continuous("Count")
    p
  })
  dev.off()
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Zheng et al.-style Plot
# Uses total plasma and PML concentrations instead of unbound/unionised
# Run a trial simulation
  mrgsimdf <- mod %>% 
    ev(time = c(0, 24, 48, 72, 96), amt = c(500, rep(250, times = 4))) %>% 
    mrgsim(end = 220, delta = 0.5) %>%
    as_tibble() %>%
    select(time, CPLAST, CMUSC, CSUBC, CPMLT) %>%
    pivot_longer(cols = contains("C"), names_to = "DVID", values_to = "DV") %>%
    mutate(DVIDf = factor(DVID, labels = c("Muscle (unb.)", "Plasma (tot.)", 
      "PML (tot.)", "Subcutis (unb.)")))
  
# Plot results of trial simulation
  pdf("evaluation/02_NONMEMValidation/GautretRegimen_ZhengPlot.pdf", width = 8, height = 10)
  p <- NULL
  p <- ggplot(mrgsimdf)
  p <- p + geom_line(aes(x = time, y = DV, colour = DVIDf), size = 1)
  p <- p + scale_x_continuous("Time (hour)", breaks = 0:9*24)
  p <- p + scale_y_log10("Concentration (ng/mL)", 
    breaks = c(1,  10, 100, 1000, 10000, 100000), labels = comma)
  p <- p + coord_cartesian(xlim = c(0, 216), ylim = c(1, 100000), expand = FALSE)
  p <- p + scale_colour_manual("Tissue", 
    values = c("darkgreen", "goldenrod", "blue", "red"))
  p <- p + theme(legend.position = "bottom")
  p 
  dev.off()
  
