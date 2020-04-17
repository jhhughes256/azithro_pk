# Typical value (PRED) plots
# -----------------------------------------------------------------------------

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Define script specific functions
# Select output function
  capture_fn <- function(x) {
    out <- c("time", "CPLAST")
    if ("CWBCT" %in% attr(mod, "capture")) { out <- c(out, c("CWBCT")) }
    else { out <- c(out, c("CPMLT")) }
    if ("AM" %in% x) {
      if ("CALMAT" %in% attr(mod, "capture")) { out <- c(out, c("CALMAT")) }
      else if ("CWBCLT" %in% attr(mod, "capture")) { out <- c(out, c("CWBCLT")) }
      else { out <- c(out, c("CPMLLT")) }
    }
    if ("Muscle" %in% x) { out <- c(out, c("CMUSC", "CSUBC")) }
    if ("Lung" %in% x) { out <- c(out, c("CLUNG")) }
    return(out)
  }
  
# Labels functions
  label_fn <- function(x) {
    out <- c("Plasma (tot.)")
    if ("CALMAT" %in% x) { out <- c(out, c("AM")) }
    if ("CLUNG" %in% x) { out <- c(out, c("Lung")) }
    if ("CMUSC" %in% x) { out <- c(out, c("Muscle", "Subcutis")) }
    if ("CPMLLT" %in% x) { out <- c(out, c("PML (ion.)")) }
    if ("CPMLT" %in% x) { out <- c(out, c("PML (tot.)")) }
    if ("CWBCLT" %in% x) { out <- c(out, c("WBC (ion.)")) }
    if ("CWBCT" %in% x) { out <- c(out, c("WBC (tot.)")) }
    out[order(out)]
  }
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Define unique regimens
  unq_regimen <- unique(digpkdf$Regimen)
  
# For each unique regimen...
  p_all <- map(unq_regimen, function(reg) {
  # Create observations dataset based on regimen
    obsdf <- digpkdf %>%
      filter(Regimen == reg) %>%
      mutate(Sample = factor(Sample)) %>% 
      mutate(Regimen = str_replace_all(Regimen, "_", " "))
  # Determine dose quantity used in regimen
    dosequant <- reg %>%
      str_extract("\\d+mg$") %>%
      str_extract("\\d+") %>%
      as.double()
  # Determine final time for simulation
    simtime <- obsdf %>%
      summarise(max(Time)) %>%
      unlist()
  # If single dose
    if (str_detect(reg, "^D1_\\d+mg$")) {
    # Simulate typical values
      preddf <- mod %>% 
        ev(amt = dosequant) %>% 
        carry_out(amt) %>%
        mrgsim_df(end = simtime + 10, delta = 0.5) %>%
        as_tibble()
  # If initial dose is given over multiple days
    } else if (str_detect(reg, "^D1-D\\d+_\\d+mg$")) {
    # Determine day of last dose (only)
      lastdose <- reg %>%
        str_extract("-D\\d+") %>%
        str_extract("\\d+") %>%
        as.double()
    # Simulate typical values
      preddf <- mod %>% 
        ev(amt = dosequant, ii = 24, addl = lastdose - 1) %>% 
        carry_out(amt) %>%
        mrgsim_df(end = simtime + 10, delta = 0.5)
  # If regimen is 500 mg D1 followed by 250 mg D2-D5
    } else {
    # Simulate typical values
      preddf <- mod %>% 
        ev(amt = c(500, rep(250, times = 4)), time = 24*0:4) %>% 
        carry_out(amt) %>%
        mrgsim_df(end = simtime + 10, delta = 0.5)
    }
  # Determine simulation outputs to capture
    sampout <- obsdf %>%
      pull(Sample) %>%
      unique() %>%
      capture_fn() %>%
      {.[. %in% names(preddf)]}
  # Subset and restructure simulation dataset
    preddf %<>%
      select_at(sampout) %>%
      pivot_longer(cols = contains("C"), names_to = "DVID", values_to = "DV") %>%
      mutate(Sample = factor(DVID, labels = label_fn(sampout)))
  # Define plot breaks
    plot_breaks <- case_when(
      simtime > 300 ~ 48,
      simtime > 100 & simtime <= 300 ~ 24,
      simtime <= 100 ~ 12)
  # Define plot limits
    plot_lims <- c(obsdf$DV, preddf$DV) %>%
      magrittr::extract(. > 0) %>%
      {c(min(.), max(.))} %>% log10() %>% magrittr::multiply_by(2) %>% 
      {c(floor(.[[1]]), ceiling(.[[2]]))} %>% magrittr::divide_by(2) %>% {10^.}
  # Define unique sample types in observed and simulated data
    sampunq <- c(levels(obsdf$Sample), levels(preddf$Sample)) %>%
      unique() %>%
      {.[order(.)]}
  # Create plot
    p <- NULL
    p <- ggplot(preddf)
    p <- p + geom_line(aes(x = time, y = DV, colour = Sample), 
      size = 1)
    p <- p + geom_errorbar(aes(x = Time, ymin = DV - SD, ymax = DV + SD, 
      colour = Sample), size = 0.5, width = simtime/50, data = obsdf)
    p <- p + geom_point(aes(x = Time, y = DV, colour = Sample, shape = Source), 
      size = 2, stroke = 1, data = obsdf)
    p <- p + scale_x_continuous("Time After First Dose (hour)", 
      breaks = 0:ceiling(simtime/plot_breaks)*plot_breaks)
    p <- p + scale_y_log10("Concentration (ng/mL, ng/mg)", 
      breaks = c(1,  10, 100, 1000, 10000, 100000, 1000000, 10000000), 
      labels = comma)
    p <- p + coord_cartesian(xlim = c(0, simtime), ylim = plot_lims)
    p <- p + scale_colour_manual("Tissue", 
      values = unlist(cPal[sampunq]))
    p <- p + scale_shape_manual("Source", 
      values = c(1, 2, 0, 5))
    p <- p + guides(colour = guide_legend(order = 1), 
      shape = guide_legend(order = 2))
    p <- p + theme(legend.position = "bottom", legend.box = "vertical", 
      legend.margin = margin())
  # Write plot to .png
    if (save_png) {
      png(paste(model_name, reg, "predplots.png", sep = "_"), width = def.width, 
        height = def.height, units = def.unit, res = def.res)
      print(p)
      dev.off()
    }
  # Add in regimen name and return plot
    p <- p + facet_wrap(~Regimen)
    p
  })
  
# Write plots to .pdf
  pdf(paste(model_name, "predplots.pdf", sep = "_"), width = def.width, 
      height = def.height)
  print(p_all)
  dev.off()
  