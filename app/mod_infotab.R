# Module UI
  
#' @title   mod_infotab_ui and mod_infotab_server
#' @description  A default shiny module
#'
#' @param id shiny id
#' @param input internal
#' @param output internal
#' @param session internal
#'
#' @rdname mod_infotab
#'
#' @keywords internal
#' @export 
#' @import shiny

mod_infotab_ui <- function(id) {
# Define namespace function for IDs
  ns <- NS(id)
# Create tagList to be used in the UI
  tagList(
    pre(htmlOutput(ns("infotext")))
  )  # tagList
}  # mod_infotab_ui

# Module Server
    
#' @rdname mod_infotab
#' @export
#' @keywords internal

mod_infotab_server <- function(input, output, session, info) {
# Define namespace function for IDs
  ns <- session$ns

# Provide session info as output
  output$infotext <- renderPrint({
    print(info)
  })  # renderPrint
  
}  # mod_infotab_server
