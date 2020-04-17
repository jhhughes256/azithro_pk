# Module UI
  
#' @title   mod_simulate_ui and mod_simulate_server
#' @description  A default shiny module
#'
#' @param id shiny id
#' @param input internal
#' @param output internal
#' @param session internal
#'
#' @rdname mod_simulate
#'
#' @keywords internal
#' @export 
#' @import shiny

mod_simulate_ui <- function(id) {
# Define namespace function for IDs
  ns <- NS(id)
# Create tagList to be used in the UI
  div(style = "display:inline-block",
    dashboardButton(ns("sim"), "Update Dosing Regimen",
      status = "success"),  # dashboardButton
    dashboardButton(ns("save"), "Save Current Output",
      status = "success")  # dashboardButton
  )  # div
}  # mod_simulate_ui

# Module Server
    
#' @rdname mod_simulate
#' @export
#' @keywords internal

mod_simulate_server <- function(input, output, session, Rinput, rv) {
# Create reactiveValues object to store simulation output
  rsim <- reactiveValues(
    out = NULL,
    save = NULL
  )
  
# Observe simulate button, when pressed:
# * save current Rinput values to rv
# * run fct_simulate_model with the input from Rinput() and store in current 
#   output `rsim$out`
  observeEvent(input$sim, {
  # Save current Rinput values to rv
    names <- c("amt", "int", "dur", "bwt", "nid")
    index <- seq(1, rv$n, by = 1)
    purrr::walk(names, ~ { rv[[.x]][index] <- na.omit(Rinput()[[.x]][index]) })
  # Simulate from model
    rsim$out <- fct_simulate_model(Rinput(), session)
  })  # observeEvent
  
# Observe save button, when pressed:
# * save current Rinput values to rv
# * save current output for comparison `rsim$save`
  observeEvent(input$save, {
  # Save current Rinput values to rv
    names <- c("amt", "int", "dur", "bwt", "nid")
    index <- seq(1, rv$n, by = 1)
    purrr::walk(names, ~ { rv[[.x]][index] <- na.omit(Rinput()[[.x]][index]) })
  # Save current simulation output to rsim
    rsim$save <- rsim$out
  })  # observeEvent
  
# Return rout object to parent module (mod_regimen_server)
  return(rsim)

}  # mod_simulate_server
