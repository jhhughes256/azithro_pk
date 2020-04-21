# Module UI
  
#' @title   mod_tablesim_ui and mod_tablesim_server
#' @description  A default Shiny module
#' Child module of mod_plotsim.
#'
#' @param id shiny id
#' @param input internal
#' @param output internal
#' @param session internal
#'
#' @rdname mod_tablesim
#'
#' @keywords internal
#' @export 
#' @import shiny

mod_tablesim_ui <- function(id) {
# Define namespace function for IDs
  ns <- NS(id)
# Create tagList to be used in the UI
  tagList(
    div(class = "col-sm-12 col-md-6", align = "center",
      tableOutput(ns("tabout"))
    ),  # div
    div(class = "col-sm-12 col-md-6", align = "center",
      tableOutput(ns("tabsave"))
    )  # div
  )  # tagList
}  # mod_tablesim_ui

# Module Server
    
#' @rdname mod_tablesim
#' @export
#' @keywords internal

mod_tablesim_server <- function(input, output, session, rsim) {
# Define namespace function for IDs
  ns <- session$ns
  
# Define function for generating application tables
  fct_tabgen <- function(rsim, label, mod) {
    rsim %>%
      purrr::when(
        !is.null(.) ~ rsim %>%
          dplyr::filter(time == max(time)) %>%
          dplyr::mutate_at(dplyr::vars(tidyselect::contains("TEC")), 
            magrittr::divide_by, 24) %>%
          dplyr::mutate_at(dplyr::vars(tidyselect::contains("TEC")), round, 1) %>%
          dplyr::mutate(AZEC50 = mod$AZEC50, AZEC90 = mod$AZEC90) %>%
          dplyr::select(
            "EC50 (ng/mL)" = AZEC50,
            "Time above EC50 in AM (days)" = ALMATEC50,
            "Time above EC50 in lung (days)" = LUNGTEC50,
            "EC90 (ng/mL)" = AZEC90,
            "Time above EC90 in AM (days)" = ALMATEC90,
            "Time above EC90 in lung (days)" = LUNGTEC90) %>%
          dplyr::mutate_all(as.character),
        is.null(.) ~ tibble::tibble(
          "EC50 (ng/mL)" = as.character(mod$AZEC50),
          "Time above EC50 in AM (days)" = "-",
          "Time above EC50 in lung (days)" = "-",
          "EC90 (ng/mL)" = as.character(mod$AZEC90),
          "Time above EC90 in AM (days)" = "-",
          "Time above EC90 in lung (days)" = "-")) %>%
      tidyr::pivot_longer(cols = tidyselect::everything(), 
        names_to = label, values_to = "Value")
  }

# Render tables to be displayed in the application
  output$tabout <- renderTable({
    fct_tabgen(rsim$out, "Current Regimen", session$userData$mod)
  })
  output$tabsave <- renderTable({
    fct_tabgen(rsim$save, "Saved Regimen", session$userData$mod)
  })

}  # mod_tablesim_server
