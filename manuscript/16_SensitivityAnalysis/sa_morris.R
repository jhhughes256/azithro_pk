# Command script for global sensitivity analysis
# ------------------------------------------------------------------------------
# Based on the Morris OAT method, where all parameters are changed from one
# random permutation to another, one at a time. This is replicated numerous
# times to determine the mean sensitivity of outcomes to the perturbation of
# a parameter when other parameters are perturbed across a range of values.
# Can only work for monotonic outcomes: AUC and Cmax qualify
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
  sa_par <- c("POPKA", "POPCLF", "POPVCF", 
    "POPQP1F", "POPVP1F", "POPQP2F", "POPVP2F", 
    "POPKIN", "POPKOUT", "POPKON", "POPKOFF", 
    "POPDFALMA", "POPDFLUNG", "POPDFMUSC", "POPDFSUBC", "POPDFWBCC")
  
# Define primary parameters
# * sa_parval - values of parameters
# * sa_k - number of parameters
# * sa_relbounds - maximum relative change of a parameter
# * sa_p - number of levels to change with the relative bounds
#     + sets the levels of the parameters used in the sensitivity analysis to:
#       p == 4: 25th percentiles; p == 6: 17th percentiles; p == 8: 12.5th 
#       percentiles of the uniform distribution with boundaries of +/- sa_relbounds
# * sa_r - number of replicates (number required depends on p)
  sa_parval <- unlist(param(mod)[sa_par])
  sa_k <- length(sa_par)
  sa_relbounds <- 0.1
  sa_p <- 8
  sa_r <- 30
  
# Define secondary parameters
# * sa_delta - standard difference between parameters
# * sa_Xi - vector of lower percentiles to be included in sensitivity analysis
#   + upper percentile is calculated as sa_Xi + sa_delta
  sa_delta <- 1/sa_p
  sa_Xi <- (0:(sa_p - 1))/(sa_p)
  
# Parameter sets for the sensitivity analysis are randomly determined using 
# matrices.
  
# B = (J_{k1,1}*X + (delta/2)*((2*B - J_{k1,k})*D + J_{k1,k}))*P
  
# Purpose of each matrix is:
# * J_{k1,1} - k+1-by-1 matrix of ones 
# * X - k samples with replacement from Xi
# * J_{k1,1} * X - k+1-by-k matrix containing the lower percentiles to used in
#     the sensitivity analysis
# * B - k+1-by-k strictly lower triangular matrix of ones (no diagonal values)
# * J_{k1,k} - k+1-by-k matrix of ones
# * 2*B - J_{k1,k} - matrix with values 1 in lower triangular matrix 
#     and -1 in diagonal and upper triangular matrix
# * D - k+1-by-k diagonal matrix with randomly samples 1 and -1 values
# * (2*B - J_{k1,k})*D - matrix with value of 1 or -1 in lower triangular matrix
#     and an inverse number in (-1 or 1) in the diagonal and upper triangular 
#     matrix
# * (sa_delta/2)*((2*B - J_{k1,k})*D + J_{k1,k}) - changes 1 to delta and -1 to 0
#     + this whole process sets up whether a parameter increases or decreases
#       by delta between simulations
# * J_{k1,1}*X + (sa_delta/2)*((2*B - J_{k1,k})*D + J_{k1,k}) - describes the
#     starting and ending percentiles for each parameter
# * P - k+1-by-k random permutation matrix (a single 1 in each row and column)
#     + Multiplication with the other matrices provides a random order for when
#       each parameter is perturbed (instead of happening in order from first
#       to final parameter)
  
# Define fixed matrices for calculating parameter sets
  sa_J_k1.1 <- matrix(1, nrow = sa_k + 1, ncol = 1)
  sa_J_k1.k <- matrix(1, nrow = sa_k + 1, ncol = sa_k)
  sa_2B <- 2*lower.tri(sa_J_k1.k, diag = FALSE)
  
# For each replicate
  sa_results <- map_dfr(seq(1, sa_r, by = 1), function(i) {
  # Define random matrices (done r times)
    sa_X <- sample(sa_Xi, sa_k, replace = TRUE)
    sa_D <- diag(2*rbinom(sa_k, 1, 0.5) - 1)
    sa_P <- diag(rep(1, sa_k))[sample(1:sa_k, sa_k, replace = FALSE), ]
    
  # Define parameter sets using defined matrices
    sa_B <- (sa_J_k1.1 %*% sa_X + (sa_delta/2)*((sa_2B - sa_J_k1.k) %*% sa_D + sa_J_k1.k)) %*% sa_P
  # For each set of parameters
    sa_simout <- map_dfr(seq(1, nrow(sa_B), by = 1), function(set) {
    # Calculate the values of the parameters to simulate from
      sa_simpar <- sa_parval + sa_parval*(2*sa_B[set,]*sa_relbounds - sa_relbounds)
    # Simulate output and capture AUC and Cmax after 1 and 21 doses
      mod %>%
        param(sa_simpar) %>%
        ev(amt = 500, ii = 24, addl = 21) %>%
        mrgsim_df(end = 528, delta = 24, Request = c("PLASAUC", "MUSCAUC", 
          "SUBCAUC", "LUNGAUC", "WBCCAUC", "ALMAAUC", "PLASCMAX", "MUSCCMAX", 
          "SUBCCMAX", "LUNGCMAX", "WBCCCMAX", "ALMACMAX")) %>%
        mutate(ID = set) %>%
        mutate_at(c("PLASAUC", "MUSCAUC", "SUBCAUC", "LUNGAUC", "WBCCAUC", 
          "ALMAAUC"), ~ c(0, diff(.x))) %>%
        filter(time %in% c(24, 528)) %>%
        mutate(time = ifelse(time == 24, "SD", "SS")) %>%
        pivot_longer(matches("^PLAS|^MUSC|^SUBC|^LUNG|^WBCC|^ALMA")) %>%
        separate(name, c("tissue", "metric"), sep = 4)
    })
  # Determine the effect of the relative change of the jth parameter for the ith replicate
  # EE_{i,j} = (Y(X_1, ..., X_i + delta_i,..., X_k) - Y(X_1, ..., X_i,..., X_k))/delta_i
  # Define delta_i
    sa_delta_i <- sa_B[sa_k + 1, ] - sa_B[1, ]
  # Determine when each parameter changes (relative change is determined from 
  #   adjacent simulations)
    sa_adj <- apply(sa_P, MARGIN = 2, function(x) {which(x == 1)})
    map_dfr(seq(1, sa_k, by = 1), function (j) {
      init_par <- sa_parval[[j]]*(1 + (sa_B[1, j]*sa_relbounds*2 - sa_relbounds))
      delta_par <- 2*sa_delta_i[[j]]*sa_relbounds*sa_parval[[j]]
      sa_simout %>%
        filter(ID %in% (sa_adj[[j]] + c(0, 1))) %>%
        mutate(par = sa_par[[j]]) %>%
        group_by(par, tissue, metric, time) %>%
        # summarise(value = diff(value)/sa_delta_i[[j]])
        summarise(Sij = diff(value)/delta_par * init_par/value[[1]])
    })
  })

  sa_summary <- sa_results %>%
    group_by(par, tissue, metric, time) %>%
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
    select(par, tissue, matches("mean$"), matches("sd$")) %>%
    write_csv("global_sa_relative.csv")
    
  p_all <- map(unique(sa_summary$tissue), function(tis) {
    p <- sa_summary %>%
      ungroup() %>%
      filter(tissue == tis & metric %in% c("AUC_SD", "AUC_SS")) %>%
      mutate(par = factor(par) %>% 
        fct_relevel(rev(c("POPKA", "POPCLF", "POPVCF", "POPQP1F", "POPVP1F", 
          "POPQP2F", "POPVP2F", "POPKIN", "POPKOUT", "POPKON", "POPKOFF", 
          "POPDFALMA", "POPDFLUNG", "POPDFMUSC", "POPDFSUBC", "POPDFWBCC")))
      ) %>% ggplot()
    p <- p + ggtitle(tis)
    p <- p + geom_point(aes(x = mSij, y = par, colour = par))
    p <- p + geom_errorbarh(aes(xmin = mSij - sdSij, xmax = mSij + sdSij, y = par, colour = par), height = 0.5)
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
  png("Global_Sensitivity_Analysis.png", width = 8, height = 10, res = 300, units = "in")
  print(p_out)
  dev.off()
  