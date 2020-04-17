# Command script for running model scripts
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
  
# Determine whether .png files have been requested
  save_png <- as.logical(read_lines("save_png.txt"))

# Load model
  model_path <- dir()[[which(str_detect(dir(), "^mod\\d+.R$"))]]
  model_name <- str_extract(model_path, "^mod\\d+")
  source(model_path)  # `mod`

# Set plotting options
# * Default plot dimensions
  def.width <- 8
  def.height <- 6
  def.unit <- "in"
  def.res <- 300

# * `ggplot2` theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.margin = unit(c(1, 1.5, 1, 1), "lines"))
  
# * Define colour palette
  cPal <- tibble(
    "Plasma (tot.)" = "#0093D0",  # pfizerblue
    "MNL" = "#A6CEE3",  # light blue
    "Muscle" = "#B2DF8A",  # light green
    "Lung" = "#33A02C",  # green
    "PML" = "#FB9A99",  # pink
    "WBC (tot.)" = "#E31A1C",  # red
    "PML (tot.)" = "#E31A1C",  # red
    "Subcutis" = "#FF7F00",  # orange
    "PBM" = "#FDBF6F",  # light orange
    "AM" = "#6A3D9A",  # purple
    "WBC (ion.)" = "#CAB2D6",  # light purple
    "PML (ion.)" = "#CAB2D6",  # light purple
    yellow = "#FFFF99",  # yellow
    brown = "#B15928"  # brown
  )
  
# Define functions
# Confidence intervals
  CI90lo <- function(x) quantile(x, probs = 0.05)
  CI90hi <- function(x) quantile(x, probs = 0.95)
  CI95lo <- function(x) quantile(x, probs = 0.025)
  CI95hi <- function(x) quantile(x, probs = 0.975)
  
# Summary functions
  vpcSummary_lst <- list(
    medianS = ~ median(.x),
    loCI90S = ~ CI90lo(.x),
    hiCI90S = ~ CI90hi(.x)
  )
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Read in external data
# * Digitised data
  digpkdf <- read_csv("AZDigitised_PK_10APR2020.csv", 
    col_types = cols(
      .default = col_character(),
      Dose = col_double(),
      Time = col_double(),
      DV = col_double(),
      SD = col_double()
  )) %>% filter(Form != "SR")
  
# * PopPK dataset used in Zhao et al. model
  poppkdf <- read_csv("A0661186_POPPK_10NOV2011.csv", na = ".",
    col_types = cols(
      .default = col_double(),
      C = col_character(),
      PROJ = col_character(),
      DATE = col_date(format = "%m/%d/%Y"),
      TIME = col_time(),
      TRTG = col_character()
  ))
   
# * PopPK + WBC dataset derived from data used in Liu et al. model
  wbcpkdf <- read_csv("A0661112_PKWBC_09APR2020.csv",
    col_types = cols(
      .default = col_double(),
      PROJ = col_character(),
      TRT = col_character(),
      TRTG = col_character(),
      UNITS = col_character(),
      C = col_character()
  ))
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Run external scripts
# Model summary
  # source("model_summary.R")
  
# Diagnostic plots
  # source("diag_plots.R")
  
# Typical value plots
  source("pred_plots.R")
  
# Model predictive check
  source("process_vpc.R")
  