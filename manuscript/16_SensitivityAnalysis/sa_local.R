# Command script for local sensitivity analysis
# ------------------------------------------------------------------------------

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare work environment
# Remove all current objects in the workspace and all plots from directory
  rm(list = ls(all = TRUE))
  unlink("*.png")
  unlink("*.pdf")
  
# Set random number generation seed
  set.seed(2684319)

# Load necessary packages for data manipulation and plotting
  library(tidyverse)
  library(mrgsolve)
  library(scales)

# Load model
  source("mod25_sa.R")  # `mod`
  
# Set plotting options
# * Default plot dimensions
  def.width <- 8
  def.height <- 6
  def.unit <- "in"
  def.res <- 300

# * `ggplot2` theme
  theme_bw2 <- theme_set(theme_bw(base_size = 10))
  theme_update(plot.margin = unit(c(1, 1.5, 1, 1), "lines"))
  
  red2green <- seq_gradient_pal(low = "#E31A1C", high = "#33A02C")
  green2blue <- seq_gradient_pal(low = "#33A02C", high = "#0093D0")
  cbPal <- c(red2green(seq(0, 1, length.out = 9)[1:8]), 
    green2blue(seq(0, 1, length.out = 8)))
  
# -----------------------------------------------------------------------------
# Morris Global Sensitivity Analysis
# Define parameters for which to calculate sensitivity
  sa_par <- c("POPKA", "POPCLF", "POPVCF", "POPQP1F", "POPVP1F", "POPQP2F", "POPVP2F", 
    "POPKIN", "POPKOUT", "POPKON", "POPKOFF", 
    "POPDFALMA", "POPDFLUNG", "POPDFMUSC", "POPDFSUBC", "POPDFWBCC")
  
# Relative changes for each parameter are taken from a uniform distribution with
#   boundaries equal to the parameter value +/- the relative change of the 
#   parameter multiplied by the parameter value.
# Define primary parameters
# * sa_parval - values of parameters
# * sa_k - number of parameters
# * sa_relbounds - maximum relative change of a parameter
# * sa_p - number of levels to change within the relative bounds
#    + p = 4 : 25th percentiles (of the uniform distribution)
#    + p = 6 : 17th percentiles
#    + p = 8 : 12.5th percentiles
#    + p = 20 : 5th percentiles
  sa_parval <- unlist(param(mod)[sa_par])
  sa_k <- length(sa_par)
  sa_relbounds <- 0.2
  sa_p <- 20
  
# Define secondary parameters
# * sa_Xi - vector of percentiles to be included in sensitivity analysis
  sa_Xi <- (0:sa_p)/(sa_p)
  print(sa_Xi*sa_relbounds*2 - sa_relbounds)
  
# For each parameter
  sa_results <- map_dfr(sa_par, function(j) {
    sa_par_j <- sa_parval[[j]]
    sa_dpar <- sa_par_j*sa_relbounds
    sa_simpar <- qunif(sa_Xi, min = sa_par_j - sa_dpar, max = sa_par_j + sa_dpar)
    sa_simout <- map_dfr(sa_simpar, function(i) {
      list(i) %>%
        set_names(j) %>%
        param(.x = mod) %>%
        ev(amt = 500, ii = 24, addl = 21) %>%
        mrgsim_df(end = 528, delta = 24, Request = c("PLASAUC", "MUSCAUC", 
          "SUBCAUC", "LUNGAUC", "WBCCAUC", "ALMAAUC", "PLASCMAX", "MUSCCMAX", 
          "SUBCCMAX", "LUNGCMAX", "WBCCCMAX", "ALMACMAX")) %>%
        mutate(ID = j, par = i) %>%
        mutate_at(c("PLASAUC", "MUSCAUC", "SUBCAUC", "LUNGAUC", "WBCCAUC", 
          "ALMAAUC"), ~ c(0, diff(.x))) %>%
        filter(time %in% c(24, 528)) %>%
        mutate(time = ifelse(time == 24, "SD", "SS")) %>%
        pivot_longer(matches("^PLAS|^MUSC|^SUBC|^LUNG|^WBCC|^ALMA")) %>%
        separate(name, c("tissue", "metric"), sep = 4)
    })
    sa_baseout <- sa_simout %>%
      filter(par == sa_par_j)
    sa_simout %>%
      filter(par != sa_par_j) %>%
      mutate(basepar = sa_par_j, basevalue = rep(sa_baseout$value, sa_p)) %>%
      mutate(Sij = (value - basevalue)/(par - basepar) * basepar/basevalue)
  })

# Create summary to determine sensitivity index
  sa_summary <- sa_results %>%
    group_by(ID, tissue, metric, time) %>%
    summarise(mSij = mean(Sij), sdSij = sd(Sij)) %>%
    arrange(tissue, metric, time, desc(abs(mSij))) %>%
    unite("metric", metric, time)
  sa_summary  %>%
    unite("value", mSij, sdSij) %>%
    pivot_wider(names_from = "metric", values_from = "value") %>%
    separate(AUC_SD, c("AUC_SD_mean", "AUC_SD_sd"), sep = "_") %>%
    separate(AUC_SS, c("AUC_SS_mean", "AUC_SS_sd"), sep = "_") %>%
    separate(CMAX_SD, c("CMAX_SD_mean", "CMAX_SD_sd"), sep = "_") %>%
    separate(CMAX_SS, c("CMAX_SS_mean", "CMAX_SS_sd"), sep = "_") %>%
    mutate_at(which(str_detect(names(.), "SD|SS")), as.double) %>%
    select(ID, tissue, matches("mean$"), matches("sd$")) %>%
    write_csv("local_sa_relative.csv")
    
  p_all <- map(unique(sa_summary$tissue), function(tis) {
    p <- sa_summary %>%
      ungroup() %>%
      filter(tissue == tis & metric %in% c("AUC_SD", "AUC_SS")) %>%
      mutate(ID = factor(ID) %>% 
        fct_relevel(rev(c("POPKA", "POPCLF", "POPVCF", "POPQP1F", "POPVP1F", 
          "POPQP2F", "POPVP2F", "POPKIN", "POPKOUT", "POPKON", "POPKOFF", 
          "POPDFALMA", "POPDFLUNG", "POPDFMUSC", "POPDFSUBC", "POPDFWBCC")))
      ) %>% ggplot()
    p <- p + ggtitle(tis)
    p <- p + geom_point(aes(x = mSij, y = ID, colour = ID))
    p <- p + geom_errorbarh(aes(xmin = mSij - sdSij, xmax = mSij + sdSij, y = ID, colour = ID), height = 0.5)
    p <- p + geom_vline(xintercept = 0, linetype = "dashed")
    p <- p + scale_colour_manual(values = cbPal)
    p <- p + theme(legend.position = "none")
    p <- p + scale_x_continuous("Sensitivity Index", limit = 1.5*c(-1, 1), 
      breaks = (-3:3)/2)
    p <- p + labs(y = NULL)
    p <- p + facet_wrap(~metric)
    p
  })
  p_out <- cowplot::plot_grid(plotlist = p_all, align = "hv", ncol = 2)
  
# Write to .png file
  png("Local_Sensitivity_Analysis.png", width = 8, height = 10, res = 300, units = "in")
  print(p_out)
  dev.off()
  