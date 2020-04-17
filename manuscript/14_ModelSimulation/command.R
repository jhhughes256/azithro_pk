# Command script for running model simulations
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
  model_path <- dir()[[which(str_detect(dir(), "^mod\\d+.R$"))]]
  model_name <- str_extract(model_path, "^mod\\d+")
  source(model_path)  # `mod`
  
# Extract omega matrix from mrgsolve model
  mrgomega <- omat(mod, make = TRUE)
	
# Define number of individuals to simulate
  nid <- 1000
	
# Define IC90 (ng/mL)
  covid_ic90 <- 6500

# Set plotting options
# * Default plot dimensions
  def.width <- 8
  def.height <- 6
  def.unit <- "in"
  def.res <- 300

# * `ggplot2` theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.margin = unit(c(1, 1.5, 1, 1), "lines"))
  
# * Determine whether .png files have been requested
  save_png <- as.logical(read_lines("save_png.txt"))
  
# Define summary functions 
# Confidence intervals
  CI90lo <- function(x) quantile(x, probs = 0.05)
  CI90hi <- function(x) quantile(x, probs = 0.95)
  CI95lo <- function(x) quantile(x, probs = 0.025)
  CI95hi <- function(x) quantile(x, probs = 0.975)
  
# Summary functions
  summary_lst <- list(
    medianS = ~ median(.x),
    loCI90S = ~ CI90lo(.x),
    hiCI90S = ~ CI90hi(.x)
  )
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Simulate regimens used in new clinical studies
# Define NCT and regimens by defining arguments for the `mrgsolve::ev` function
  regimen_list <- list(
    list(id = "D1_1000mg",
      label = "Day 1: 1000 mg",
      nct = c("NCT04332107"),
      reg = list(amt = 1000) ),
    list(id = "D1-3_500mg",
      label = "Days 1-3: 500 mg",
      nct = c("NCT04332094"),
      reg = list(amt = 500, ii = 24, addl = 2) ),
    list(id = "D1-5_500mg",
      label = "Days 1-5: 500 mg",
      nct = c("NCT04329572"),
      reg = list(amt = 500, ii = 24, addl = 4) ),
    list(id = "D1-7_500mg",
      label = "Days 1-7: 500 mg",
      nct = c("NCT04322123"),
      reg = list(amt = 500, ii = 24, addl = 6) ),
    list(id = "D1-10_500mg",
      label = "Days 1-10: 500 mg",
      nct = c("NCT04321278"),
      reg = list(amt = 500, ii = 24, addl = 9) ),
    list(id = "D1_500mg_D2-5_250mg",
      label = "Day 1: 500 mg; Days 2-5: 250 mg",
      nct = c("NCT04335552", "NCT04336332", "NCT04329832", "NCT04334382", "NCT04338698", "NCT04324463"),
      reg = list(amt = c(500, rep(250, times = 4)),
      time = c(seq(0, 24*4, by = 24))) ),
    list(id = "D1_500mg_D2-7_250mg",
      label = "Day 1: 500 mg; Days 2-7: 250 mg",
      nct = c("NCT04328272"),
      reg = list(amt = c(500, rep(250, times = 6)),
      time = c(seq(0, 24*6, by = 24))) )
  )
  
# Simulate ETA values and weights for `nid` subjects
  mrgidata <- mrgomega %>%
    {MASS::mvrnorm(n = nid, mu = rep(0, times = nrow(.)), Sigma = .)} %>%
    as_tibble() %>%
    set_names(paste0("ETA", 1:nrow(mrgomega))) %>%
    add_column(ID = 1:(nid), .before = "ETA1")
  
# For each regimen
  p_all <- map(regimen_list, function(regimen) {
    out <- regimen %>%
      pluck("reg") %>%
      list_modify(x = mod) %>%
      {do.call("ev", .)} %>%
      idata_set(mrgidata) %>%
      carry_out(amt) %>%
      mrgsim_df(end = 490, delta = 0.5) %>%
      select(ID, time, CALMAT, CLUNG, CPLAST, CWBCT) %>%
      pivot_longer(cols = contains("C"), names_to = "DVID", values_to = "DV") %>%
      mutate(Sample = factor(DVID, labels = c("AM", "Lung", "Plasma", "WBC"))) %>%
      group_by(Sample, time) %>%
      summarise_at("DV", summary_lst) %>%
      mutate(Regimen = pluck(regimen, "label"))
      
    p <- ggplot(out)
    p <- p + geom_line(aes(x = time/24, y = medianS, colour = Sample), size = 1)
    p <- p + geom_ribbon(aes(x = time/24, ymin = loCI90S, ymax = hiCI90S, fill = Sample), alpha = 0.3)
    p <- p + geom_hline(yintercept = covid_ic90, linetype = "dashed")
    p <- p + scale_x_continuous("\nTime after First Dose (days)",
      breaks = seq(0, 480/24, by = 48/24))
    p <- p + scale_y_log10("Concentration (ng/mL)", 
      breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000), labels = comma)
    p <- p + scale_colour_manual("Tissue", 
      values = c("#6A3D9A", "#33A02C", "#0093D0", "#E31A1C"))
    p <- p + scale_fill_manual("", guide = FALSE,
      values = c("#6A3D9A", "#33A02C", "#0093D0", "#E31A1C"))
    p <- p + coord_cartesian(xlim = c(10/24, 480/24), ylim = c(1, 3000000))
    p <- p + theme(legend.position = "bottom")
    if (save_png) {
      png(paste(model_name, pluck(regimen, "id"), "AZsim.png", sep = "_"), width = def.width, 
        height = def.height, units = def.unit, res = def.res)
      print(p)
      dev.off()
    }
    return(p + facet_wrap(~Regimen))
  })
  
  pdf(paste(model_name, "AZsim.pdf", sep = "_"), width = def.width, height = def.height)
  print(p_all)
  dev.off()

# Remove all current objects in the workspace
  rm(list = ls(all = TRUE))