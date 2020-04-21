# Module UI
  
#' @title   mod_plotsim_ui and mod_plotsim_server
#' @description  Shiny module which handles display of concentration-time plot.
#' Includes data processing required to generate the lines and ribbons in the
#' plot. Parent module for mod_tablesim.
#'
#' @param id shiny id
#' @param input internal
#' @param output internal
#' @param session internal
#'
#' @rdname mod_plotsim
#'
#' @keywords internal
#' @export 
#' @import shiny

mod_plotsim_ui <- function(id) {
# Define namespace function for IDs
  ns <- NS(id)
# Create tagList to be used in the UI
# * box width according to parent div
  box(width = "100%", title = "Azithromycin Concentrations", align = "center",
    plotOutput(ns("mainplot"), height = "450px"),
    footer = mod_tablesim_ui(ns("table")),
    status = "success", solidHeader = TRUE
  )  # box
}  # mod_plotsim_ui

# Module Server
    
#' @rdname mod_plotsim
#' @export
#' @keywords internal

mod_plotsim_server <- function(input, output, session, rsim) {
# Define namespace function for IDs
  ns <- session$ns
  
# Create function for generating plotted lines from simulation data
# * function is applied to both current and saved outputs
# * removes TEC50 and TEC90 columns prior to processing
# * plot_summary_fn can be found in utils.R
  fct_plotlines <- function(rsim, plot_summary_arg) {
    rsim %>%
      dplyr::select(-tidyselect::contains("TEC")) %>%
      dplyr::mutate(time = time/24) %>%
      tidyr::pivot_longer(cols = -c("ID", "time"),
        names_to = "Tissue", values_to = "DV") %>%
      dplyr::mutate(Tissue = factor(Tissue, 
        labels = c("AM", "Lung", "Plasma", "WBC"))) %>%
      dplyr::group_by(time, Tissue) %>%
      purrr::when(
        length(unique(.$ID)) == 1 ~ dplyr::summarise(., median = median(DV)),
        length(unique(.$ID)) > 1 ~ dplyr::summarise_at(., "DV", 
          plot_summary_fn(plot_summary_arg))
      ) %>%
      dplyr::ungroup()
  }
  
# Generate plotted line from current and saved simulation data
# * generate 90%, 80%, 60%, 40% and 20% prediction intervals for current output
# * generate 90% prediction intervals for saved output
  Rlines <- reactive({
    fct_plotlines(rsim$out, c(0.9, 0.8, 0.6, 0.4, 0.2))
  })
  Slines <- reactive({
    fct_plotlines(rsim$save, c(0.9))
  })
  
# Create function for extracting prediction intervals from summary data
  fct_plotribbons <- function(lines) {
    lines %>%
      dplyr::select(time, Tissue, tidyselect::contains("pi")) %>%
      tidyr::pivot_longer(cols = tidyselect::contains("pi"), 
        names_to = "PI", values_to = "value") %>%
      tidyr::separate(col = "PI", into = c("PI", "id"), sep = -2) %>%
      tidyr::pivot_wider(id_cols = c("time", "Tissue", "PI"),
        names_from = "id", values_from = "value") %>%
      dplyr::mutate(PI = paste(Tissue, PI))
  }
  
# Generate plotted ribbon from current and saved summary data
# * try used when using fct_plotribbons on Slines() to help protect from error
#     + class(Sribbon()) == "try-error" when saved output is from a simulation
#       with `rv$nid` == 1, while current output is from a simulation with
#       `rv$nid` > 1.
  Rribbon <- reactive({
    fct_plotribbons(Rlines())
  })
  Sribbon <- reactive({
    try(fct_plotribbons(Slines()))
  })

# Create main plot to be displayed in the application
  output$mainplot <- renderPlot({
  # Define azithromycin IC90 (ng/mL)
    az_ic90 <- 6500  # TODO: do we want the user to access this?
  # Only produce plot if simulation has occurred
    if (!is.null(rsim$out)) {  # if simulation has occurred
    # Define plot
      p <- ggplot2::ggplot(data = Rlines())
    # Azithromycin Model current predictions
      p <- p + ggplot2::geom_line(ggplot2::aes(x = time, y = median, 
        colour = Tissue), size = 1)
      if (length(unique(rsim$out$ID)) > 1) {  # If nid > 1
        p <- p + ggplot2::geom_ribbon(ggplot2::aes(x = time, ymin = lo, ymax = hi, 
          group = PI, fill = Tissue), alpha = 0.1, data = Rribbon())
      }
    # Azithromycin Model saved predictions
    # Display if simulation output has been saved, and isn't identical to current output
      if (!is.null(rsim$save) & !isTRUE(all.equal(rsim$save, rsim$out))) {
        p <- p + ggplot2::geom_line(ggplot2::aes(x = time, y = median, 
          colour = Tissue), size = 1, linetype = "dashed", data = Slines())
        if (length(unique(rsim$out$ID)) > 1 & !"try-error" %in% class(Sribbon())) {  # If nid > 1
          p <- p + ggplot2::geom_line(ggplot2::aes(x = time, y = lo, 
            colour = Tissue), size = 1, alpha = 0.7, linetype = "dotted", 
            data = Sribbon())
          p <- p + ggplot2::geom_line(ggplot2::aes(x = time, y = hi, 
            colour = Tissue), size = 1, alpha = 0.7, linetype = "dotted", 
            data = Sribbon())
        }
      }
    # Azithromycin SARS-CoV-2 IC90
      p <- p + ggplot2::geom_hline(yintercept = az_ic90, linetype = "dashed")
      p <- p + ggplot2::geom_text(x = 21, y = log10(az_ic90) + 0.2, label = "EC90")
    # Axes
      p <- p + ggplot2::scale_y_log10(paste("Concentration (ng/mL)"),
        labels = scales::comma)
      p <- p + ggplot2::scale_x_continuous("Time After First Dose (days)",
        breaks = seq(min(rsim$out$time), max(rsim$out$time), by = 7),
        minor_breaks = seq(min(rsim$out$time), max(rsim$out$time), by = 1))
    # Legend
      p <- p + ggplot2::scale_colour_manual("Tissue", 
        values = c("#6A3D9A", "#33A02C", "#0093D0", "#E31A1C"))
      p <- p + ggplot2::scale_fill_manual("Tissue", 
        values = c("#6A3D9A", "#33A02C", "#0093D0", "#E31A1C"))
      p <- p + ggplot2::theme(legend.position = "bottom", 
        legend.margin = ggplot2::margin())
      return(p)
    }
  })
  
# Call module responsible for table output and pass in inputs
# Assign output from module to rsim
  callModule(mod_tablesim_server, "table", rsim = rsim)
  
}  # mod_plotsim_server
