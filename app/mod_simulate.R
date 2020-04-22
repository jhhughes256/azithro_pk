# Module UI
  
#' @title   mod_simulate_ui and mod_simulate_server
#' @description  Shiny module which handles all dosing regimen input. Includes 
#' buttons to start simulation and save/clear simulation output for comparison.
#' Child module of mod_regimen.
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
# * dashboardButton wrapper function can be found in `utils.R`
  div(align = "center",
    dashboardButton(ns("sim"), "Update Dosing Regimen",
      status = "warning", width = "100%"),  # dashboardButton
    br(), br(),  # spacing
    fluidRow(
      column(6,
        dashboardButton(ns("save"), "Save Current Output",
          status = "warning", width = "100%")  # dashboardButton
      ),
      column(6,
        dashboardButton(ns("clear"), "Clear Saved Output",
          status = "warning", width = "100%")  # dashboardButton
      )
    )  # div
  )  # tagList
}  # mod_simulate_ui

# Module Server
    
#' @rdname mod_simulate
#' @export
#' @keywords internal

mod_simulate_server <- function(input, output, session, Rinput, rv) {
# Create reactiveValues object to store simulation output
# * `rsim$out` contains the current output from mrgsolve simulation
# * `rsim$save` contains the saved output from mrgsolve simulation
  rsim <- reactiveValues(
    out = NULL,
    save = NULL
  )
  
# Observe simulate button, when pressed:
# * check to make sure that Rinput() has been called successfully
# * if regimen is user-defined save current Rinput() values to `rv`
# * run fct_simulate_model with the input from Rinput() and store in current 
#   output `rsim$out`
# * if Rinput() wasn't called successfully, let user know that they clicked the
#   simulation button a bit too quickly and should try again.
  observeEvent(input$sim, {
    if (!"try-error" %in% class(Rinput())) {
    # Save current Rinput values to rv if required
      # if (Rinput()$choose == "0") {
      #   names <- c("amt", "int", "dur")
      #   index <- seq(1, rv$n, by = 1)
      #   purrr::walk(names, ~ { rv[[.x]][index] <- na.omit(Rinput()[[.x]][index]) })
      # }
    # Simulate from model
      rsim$out <- fct_simulate_model(Rinput(), session)
    } else {
      showNotification("Sorry you were a bit quick! Try again.", type = "warning")
    }
  })  # observeEvent
  
# Observe save button, when pressed save current output to `rsim$save`
  observeEvent(input$save, {
    rsim$save <- rsim$out
  })  # observeEvent
  
# Observe clear button, when pressed clear current output from `rsim$save`
  observeEvent(input$clear, {
    rsim$save <- NULL
  })  # observeEvent
  
# Return rout object to parent module (mod_regimen_server)
  return(rsim)

}  # mod_simulate_server
