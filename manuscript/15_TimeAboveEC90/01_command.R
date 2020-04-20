# Obtain numbers for manuscript
# ------------------------------------------------------------------------------
# * Determine comparative concentrations in tissues at 24 hours
# * Determine time above IC90 in AM and lung tissue for each regimen
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

# Load model
  source("mod25_alt.R")  # `mod`
  
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
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Determine comparative concentrations in tissues at 24 hours
  map(regimen_list[2], function(regimen) {
    out <- regimen %>%
      pluck("reg") %>%
      list_modify(x = mod) %>%
      {do.call("ev", .)} %>%
      carry_out(amt) %>%
      mrgsim_df(end = 48, delta = 0.5, Request = c("CWBCT", "CALMAT")) %>%
      filter(time == 24) %>%
      with(CALMAT/CWBCT)
  })
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Determine time above EC50/EC90 in AM and lung tissue for each regimen
  map(regimen_list, function(regimen) {
    out <- regimen %>%
      pluck("reg") %>%
      list_modify(x = mod) %>%
      {do.call("ev", .)} %>%
      carry_out(amt) %>%
      mrgsim_df(end = 480, delta = 0.5, 
        Request = c("ALMATEC50", "LUNGTEC50", "ALMATEC90", "LUNGTEC90")) %>%
      filter(time == 480) %>%
      mutate_if(str_detect(names(.), "TEC"), magrittr::divide_by, 24) %>%
      mutate_if(str_detect(names(.), "TEC"), round, 1)
  })
      
