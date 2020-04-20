# Module UI
  
#' @title   mod_plotsim_ui and mod_plotsim_server
#' @description  A default shiny module
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
  box(width = NULL, title = "Azithromycin Concentrations", align = "center",
    plotOutput(ns("mainplot"), height = "450px"),
    status = "primary", solidHeader = TRUE
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
  lines_fn <- function(rsim, plot_summary_arg) {
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
  Rlines <- reactive({
    lines_fn(rsim$out, c(0.9, 0.8, 0.6, 0.4, 0.2))
  })
  Slines <- reactive({
    lines_fn(rsim$save, c(0.9))
  })
  
# Create function for generating prediction intervals from simulation data
  ribbon_fn <- function(lines) {
    lines %>%
      dplyr::select(time, Tissue, tidyselect::contains("pi")) %>%
      tidyr::pivot_longer(cols = tidyselect::contains("pi"), 
        names_to = "PI", values_to = "value") %>%
      tidyr::separate(col = "PI", into = c("PI", "id"), sep = -2) %>%
      tidyr::pivot_wider(id_cols = c("time", "Tissue", "PI"),
        names_from = "id", values_from = "value") %>%
      dplyr::mutate(PI = paste(Tissue, PI))
  }
  
# Generate plotted ribbon from current and saved simulation data
  Rribbon <- reactive({
    ribbon_fn(Rlines())
  })
  Sribbon <- reactive({
    try(ribbon_fn(Slines()))
  })

# Create main 
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
      p <- p + ggplot2::geom_text(x = 21, y = log10(az_ic90) + 0.2, label = "IC90")
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
  
}  # mod_plotsim_server
