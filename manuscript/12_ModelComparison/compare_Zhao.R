# Model comparison: Zheng vs. Zhao
# -----------------------------------------------------------------------------
# Comparison of typical values for Zheng et al. and Zhao et al. (Pfizer) models,
# to highlight the magnitude of difference between the models.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Prepare work environment
# Remove all current objects in the workspace
  rm(list = ls(all = TRUE))

# Load necessary packages for data manipulation and plotting
  library(tidyverse)
  library(mrgsolve)
  library(scales)

# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.margin = unit(c(1, 1.5, 1, 1), "lines"))

# Load external files
# * Compile models
  source("model.R")
  source("modelBallow.R")
  source("modelLiu.R")
  source("modelSalman.R")
  source("modelZhang.R")
  source("modelZhao.R")
  
# - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - - 
# Compare PRED values from both models
# Run a trial simulation for Zheng et al.
  zhengdf <- mod %>% 
    ev(amt = c(500, rep(250, times = 4)), time = 24*0:4) %>% 
    mrgsim(end = 368, delta = 0.5) %>%
    as_tibble()
  
# Run a trial simulation for Ballow et al.
  balldf <- modBall %>% 
    ev(amt = c(500, rep(250, times = 4)), time = 24*0:4) %>% 
    mrgsim(end = 368, delta = 0.5) %>%
    as_tibble()
  
# Run a trial simulation for Muto et al.
  liupdf <- modLiu %>% 
    ev(amt = c(500, rep(250, times = 4)), time = 24*0:4) %>% 
    mrgsim(end = 368, delta = 0.5) %>%
    as_tibble()
  
# Run a trial simulation for Zhao et al.
  salmdf <- modSalm %>% 
    param(PREG = 0, WT = 78) %>%
    ev(amt = c(500, rep(250, times = 4)), time = 24*0:4) %>% 
    mrgsim(end = 368, delta = 0.5) %>%
    as_tibble()
  
# Run a trial simulation for Zhang et al.
  zhangdf <- modZhang %>% 
    param(AGE = 29, WT = 78) %>%
    ev(amt = c(500, rep(250, times = 4)), time = 24*0:4) %>% 
    mrgsim(end = 368, delta = 0.5) %>%
    as_tibble()
  
# Run a trial simulation for Zhao et al.
  zhaodf <- modZhao %>% 
    param(WT = 78) %>%
    ev(amt = c(500, rep(250, times = 4)), time = 24*0:4) %>% 
    mrgsim(end = 368, delta = 0.5) %>%
    as_tibble()
  
# Plot results of trial simulation
  pdf("compare_typical.pdf", width = 8, height = 6)
  p <- NULL
  p <- ggplot()
  p <- p + geom_line(aes(x = time, y = CPLAST, colour = "Zheng et al."), 
    size = 1, data = zhengdf)
  p <- p + geom_line(aes(x = time, y = CPLAST, colour = "Ballow et al."), 
    size = 1, data = balldf)
  p <- p + geom_line(aes(x = time, y = CPLAST, colour = "Liu et al."), 
    size = 1, data = liupdf)
  p <- p + geom_line(aes(x = time, y = CPLAST, colour = "Salman et al."), 
    size = 1, data = salmdf)
  p <- p + geom_line(aes(x = time, y = CPLAST, colour = "Zhang et al."), 
    size = 1, data = zhangdf)
  p <- p + geom_line(aes(x = time, y = CPLAST, colour = "Zhao et al."), 
    size = 1, data = zhaodf)
  p <- p + scale_x_continuous("Time (hour)", breaks = 0:21*24)
  p <- p + scale_y_log10("Concentration (ng/mL)", 
    breaks = c(1,  10, 100, 1000, 10000, 100000, 1000000), labels = comma)
  p <- p + coord_cartesian(xlim = c(0, 360), ylim = c(1, 1000))
  p <- p + scale_colour_manual("Model", 
    values = c("#0093D0", "#D55E00", "#009E73", "#E69F00", "#56B4E9", "#CC79A7"))
  p <- p + theme(legend.position = "bottom")
  p  
  dev.off()
  
  png("compare_typical.png", width = 12, height = 6, units = "in", res = 300)
  p
  dev.off()