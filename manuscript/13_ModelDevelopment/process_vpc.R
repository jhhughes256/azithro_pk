# Process visual predictive check plots
# -----------------------------------------------------------------------------

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Define VPC functions 
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
# Simulate A0661112 profiles 1000 times
# Define number of simulations
  nsim <- 1000
  
# Prepare observation data for mrgsolve
  mrglist <- list(
    poppkdf %>%
      filter(is.na(C) & (FLAG %in% c(0, 100) | (FLAG == 2 & AMT == DOS1))) %>%
      mutate(TAFDBIN = round(TAFD)) %>%
      select(PROT, ID, AMT, TIME = TAFD, TAD, TAFDBIN, NTPD, DV, MDV, EVID, CMT, WT = BWT),
    wbcpkdf %>% 
      filter(is.na(C) & TRT == "B" & FLAG %in% c(0, 1)) %>%
      mutate(DV = DV*1000, CMT = 1, TAFDBIN = round(TAFD)) %>%
      select(PROT, NSID, AMT, TIME = TAFD, TAFDBIN, NTPD, DV, MDV, EVID, CMT, FLAG, WT = BWT, AGE) %>%
      mutate(ID = as.double(factor(NSID)))
  )
# Extract omega matrix from mrgsolve model
	mrgomega <- omat(mod, make = TRUE)
  
# For each set of observation data
  p_all <- map(mrglist, function(mrgdata) {
  # Define number of individuals per simulation
  	nid <- length(unique(mrgdata$ID))
  # Create dataset specifically for Liu et al. VPC
  	if ("NSID" %in% names(mrgdata)) {
	  # Prepare simulation input data for mrgsolve (24-48 hours)
      impdata <- mrgdata %>%
        filter(MDV == 0) %>%
        distinct(ID, .keep_all = TRUE) %>%
        mutate(DV = 0) %>%
        nest(time_cols = c(TIME, NTPD, TAFDBIN)) %>%
        mutate(time_cols = map(time_cols, ~ tibble(
          NTPD = c(0.5, 1, 2, 3, 4, 6, 8, 12, 23.5), 
          TIME = 24 + NTPD, TAFDBIN = TIME)
        )) %>%
        unnest(time_cols)
  	}
  # Prepare individual data for mrgsolve
    mrgidata <- mrgomega %>%
      {MASS::mvrnorm(n = nid*nsim, mu = rep(0, times = nrow(.)), Sigma = .)} %>%
      as_tibble() %>%
      set_names(paste0("ETA", 1:nrow(mrgomega))) %>%
      add_column(ID = 1:(nid*nsim), .before = "ETA1")
  # Simulate profiles
    mrgdf <- mod %>%
      data_set(map_dfr(1:nsim, function(sim) {
        if ("NSID" %in% names(mrgdata)) {
          bind_rows(
            mutate(mrgdata, SIM = sim, ID = ID + nid*(sim - 1)),
            mutate(impdata, SIM = sim, ID = ID + nid*(sim - 1))
          ) %>%
          arrange(SIM, ID, TIME, desc(AMT))
        } else {
          mutate(mrgdata, SIM = sim, ID = ID + nid*(sim - 1))
        }
      })) %>%
      idata_set(mrgidata) %>%
      carry_out(SIM, NSID, AMT, DV, TAFDBIN) %>%
      mrgsim() %>%
      as_tibble()
  # Summarise total plasma concentration for each SIM and TAFDBIN
    summdf <- mrgdf %>%
      group_by(SIM, TAFDBIN) %>%
      summarise_at("CPLAST", vpcSummary_lst)
  # Create plot
    p <- NULL
    p <- ggplot(data = mrgdata)
  # Define observed elements
    p <- p + geom_point(aes(x = TAFDBIN, y = DV), colour = "#0093D0", 
      shape = 1)
    p <- p + stat_summary(aes(x = TAFDBIN, y = DV), fun.y = median,
      geom = "line", colour = "#E31A1C", size = 1)
    p <- p + stat_summary(aes(x = TAFDBIN, y = DV), fun.y = CI90lo,
      geom = "line", colour = "#E31A1C", linetype = "dashed", size = 1)
    p <- p + stat_summary(aes(x = TAFDBIN, y = DV), fun.y = CI90hi,
      geom = "line", colour = "#E31A1C", linetype = "dashed", size = 1)
  # Define simulated elements
    p <- p + stat_summary(aes(x = TAFDBIN, y = medianS), 
      fun.ymin = "CI95lo", fun.ymax = "CI95hi",
      data = summdf, geom = "ribbon", alpha = 0.3, fill = "#E31A1C")
    p <- p + stat_summary(aes(x = TAFDBIN, y = medianS), fun.y = median, 
      data = summdf, geom = "line", colour = "black", size = 1)
    p <- p + stat_summary(aes(x = TAFDBIN, y = loCI90S),
      fun.ymin = "CI95lo", fun.ymax = "CI95hi",
      data = summdf, geom = "ribbon", alpha = 0.3, fill = "#0093D0")
    p <- p + stat_summary(aes(x = TAFDBIN, y = loCI90S), fun.y = median, 
      data = summdf, geom = "line", colour = "black", 
      linetype = "dashed", size = 1)
    p <- p + stat_summary(aes(x = TAFDBIN, y = hiCI90S), 
      fun.ymin = "CI95lo", fun.ymax = "CI95hi",
      data = summdf, geom = "ribbon", alpha = 0.3, fill = "#0093D0")
    p <- p + stat_summary(aes(x = TAFDBIN, y = hiCI90S), fun.y = median, 
      data = summdf, geom = "line", colour = "black", 
      linetype = "dashed", size = 1)
  # Finalise plot
    p <- p + scale_y_log10(paste0("Plasma Concentration (ng/mL)\n"), 
      breaks = c(3, 10, 30, 100, 300, 1000, 3000))
    p <- p + scale_x_continuous(paste0("\nTime after First Dose (hours)"), 
      breaks = 0:9*24)
    p <- p + coord_cartesian(ylim = c(3, 1000))
  # Write plot to .png
    if (save_png) {
      png(paste(model_name, paste0("A066", unique(mrgdata$PROT)), "PKVPC.png", sep = "_"), 
        width = def.width, height = def.height, units = def.unit, res = def.res)
      print(p)
      dev.off()
    }
    p <- p + facet_wrap(~PROT, labeller = function(x) paste0("A066", x))
    return(p)
  })
  
  pdf(paste(model_name, "PKVPC.pdf", sep = "_"), width = def.width, height = def.height) 
  print(p_all)
  dev.off()

  if (any(c("CPMLT", "CWBCT") %in% attr(mod, "capture"))) {
  # Filter data for WBC concentrations
    mrgdata <- wbcpkdf %>% 
      filter(is.na(C) & TRT == "B" & FLAG %in% c(0, 100, 101)) %>%
      mutate(DV = DV*1000, CMT = 1, TAFDBIN = round(TAFD)) %>%
      select(PROT, NSID, AMT, TIME = TAFD, NTPD, TAFDBIN, DV, MDV, EVID, CMT, FLAG, WT = BWT, AGE) %>%
      mutate(ID = as.double(factor(NSID)))
  # Define number of individuals per simulation
    nid <- length(unique(mrgdata$ID))
  # Prepare individual data for mrgsolve
    mrgidata <- mrgomega %>%
      {MASS::mvrnorm(n = nid*nsim, mu = rep(0, times = nrow(.)), Sigma = .)} %>%
      as_tibble() %>%
      set_names(paste0("ETA", 1:nrow(mrgomega))) %>%
      add_column(ID = 1:(nid*nsim), .before = "ETA1")
  # Simulate profiles
    mrgdf <- mod %>%
      data_set(map_dfr(1:nsim, ~ mutate(mrgdata, SIM = .x, ID = ID + nid*(.x - 1)))) %>%
      idata_set(mrgidata) %>%
      carry_out(SIM, NSID, AMT, DV, NTPD, TAFDBIN, FLAG) %>%
      mrgsim_df()
  # Summarise total cell concentration for each SIM and NTPD
    cell_pred <- if_else("CWBCT" %in% attr(mod, "capture"), "CWBCT", "CPMLT")
    summdf_ntpd <- mrgdf %>%
      group_by(SIM, NTPD) %>%
      summarise_at(cell_pred, vpcSummary_lst)
    summdf_tafdbin <- mrgdf %>%
      group_by(SIM, TAFDBIN) %>%
      summarise_at(cell_pred, vpcSummary_lst)
  # Create wbc VPC plot function
    wbcvpc_fn <- function(mrgdata, summdf, time) {
    # Define x-axis label from time value
      x_label <- case_when(
        time == "NTPD" ~ "Time after Dose",
        time == "TAFDBIN" ~ "Time after First Dose")
    # Define plot observed data data
      p <- ggplot(data = filter(mrgdata, FLAG %in% c(0, 100, 101)))
    # Define observed elements
      p <- p + geom_point(aes_string(x = time, y = "DV"), colour = "#0093D0",
        shape = 1)
      p <- p + stat_summary(aes_string(x = time, y = "DV"), fun.y = median,
        geom = "line", colour = "#E31A1C", size = 1)
      p <- p + stat_summary(aes_string(x = time, y = "DV"), fun.y = CI90lo,
        geom = "line", colour = "#E31A1C", linetype = "dashed", size = 1)
      p <- p + stat_summary(aes_string(x = time, y = "DV"), fun.y = CI90hi,
        geom = "line", colour = "#E31A1C", linetype = "dashed", size = 1)
    # Define simulated elements
      p <- p + stat_summary(aes_string(x = time, y = "medianS"),
        fun.ymin = "CI95lo", fun.ymax = "CI95hi",
        data = summdf, geom = "ribbon", alpha = 0.3, fill = "#E31A1C")
      p <- p + stat_summary(aes_string(x = time, y = "medianS"), fun.y = median,
        data = summdf, geom = "line", colour = "black", size = 1)
      p <- p + stat_summary(aes_string(x = time, y = "loCI90S"),
        fun.ymin = "CI95lo", fun.ymax = "CI95hi",
        data = summdf, geom = "ribbon", alpha = 0.3, fill = "#0093D0")
      p <- p + stat_summary(aes_string(x = time, y = "loCI90S"), fun.y = median,
        data = summdf, geom = "line", colour = "black",
        linetype = "dashed", size = 1)
      p <- p + stat_summary(aes_string(x = time, y = "hiCI90S"),
        fun.ymin = "CI95lo", fun.ymax = "CI95hi",
        data = summdf, geom = "ribbon", alpha = 0.3, fill = "#0093D0")
      p <- p + stat_summary(aes_string(x = time, y = "hiCI90S"), fun.y = median,
        data = summdf, geom = "line", colour = "black",
        linetype = "dashed", size = 1)
    # Finalise plot
      p <- p + scale_y_log10(paste("WBC Concentration (ng/mL)\n"),
        breaks = c(300, 1000, 3000, 10000, 30000, 100000, 300000), labels = comma)
      p <- p + scale_x_continuous(paste0("\n", x_label, " (hours)"),
        breaks = 0:9*24)
      p <- p + coord_cartesian(ylim = c(2000, 300000))
    # Write VPC to .png file
      if (save_png) {
        png(paste(model_name, "WBCVPC", paste0(time, ".png"), sep = "_"), width = def.width,
          height = def.height, units = def.unit, res = def.res)
        print(p)
        dev.off()
      }
    # Return plot
      return(p)
    }
  # Create list of VPC plot using function
    p_wbc <- list(
      wbcvpc_fn(mrgdata, summdf_ntpd, "NTPD"), 
      wbcvpc_fn(mrgdata, summdf_tafdbin, "TAFDBIN"))
  
  # Save as .pdf file
    pdf(paste(model_name, "WBCVPC.pdf", sep = "_"), width = def.width, height = def.height) 
    print(p_wbc)
    dev.off()
  }
  
# Remove all current objects in the workspace
  rm(list = ls(all = TRUE))