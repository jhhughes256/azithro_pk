# Module UI
  
#' @title   mod_regimen_ui and mod_regimen_server
#' @description  A default shiny module
#'
#' @param id shiny id
#' @param input internal
#' @param output internal
#' @param session internal
#'
#' @rdname mod_regimen
#'
#' @keywords internal
#' @export 
#' @import shiny

mod_regimen_ui <- function(id) {
# Define namespace function for IDs
  ns <- NS(id)
# Create tagList to be used in the UI
  box(width = NULL, title = "Regimen Information", align = "center",
    numericInput(ns("bwt"), "Patient Body Weight (kg):",
      value = 79, min = 1, step = 0.1),
    selectInput(ns("choose"), "Select dosing regimen:",
      choices = c(
        "Days 1-3: 500 mg" = 1,
        "Day 1: 500mg; Days 2-5: 250 mg" = 2,
        "Day 1: 1000mg" = 3,
        "Days 1-10: 500mg" = 4,
        "Custom" = 0
      ), selected = 1),
    uiOutput(ns("regimen")),
    br(),
    div(style = "display:in-line",
      actionButton(ns("add"), "Add"),
      actionButton(ns("rem"), "Remove")
    ),  # div
    br(),
    sliderInput(ns("nid"), "Number of Individuals for Simulation:",
      value = 1, min = 1, max = 1000, step = 1, width = "90%"),
    footer = mod_simulate_ui(paste(id, "sim", sep = "-")),  # gotta rememeber to add parent to namespace for nested modules!!!
    status = "success", solidHeader = TRUE
  )  # box
}  # mod_regimen_ui

# Module Server
    
#' @rdname mod_regimen
#' @export
#' @keywords internal

mod_regimen_server <- function(input, output, session) {
# Define namespace function for IDs
  ns <- session$ns
  
# Create reactiveValues object to store input values
  rv <- reactiveValues(
    n = 1,
    amt = rep(500, 3),
    int = rep(24, 3),
    dur = rep(3, 3),
    bwt = 79,
    nid = 1
  )
  
# Define reactive function for dynamic ui
  Rui <- reactive({
    purrr::map(seq(1, rv$n, by = 1), function(i) {
      if (i == 1) { 
        ilab <- c("Dose:", "Interval:", "Duration:") 
        padding <- "padding-top:20px; padding-left:0px"
      } else { 
        ilab <- NULL 
        padding <- "padding-left:0px"
      }
      column(12, fluidRow(
        column(3,
          numericInput(ns(paste0("amt", i)), 
            label = ilab[[1]], 
            value = rv$amt[[i]])
        ), # column
        column(1,
          style = padding,
          h5(strong("mg, every"))
        ), # column
        column(3,
          numericInput(ns(paste0("int", i)), 
            label = ilab[[2]], 
            value = rv$int[[i]])
        ), # column
        column(1,
          style = padding,
          h5(strong("hours, for"))
        ), # column
        column(3,
          numericInput(ns(paste0("dur", i)), 
            label = ilab[[3]], 
            value = rv$dur[[i]])
        ), # column
        column(1,
          style = padding,
          h5(strong("days"))
        )  # column
      ))  # column.fluidRow
    })  # map
  })  # reactive
  
# Define reactive function containing input values
  Rinput <- reactive({
    list(
      amt = purrr::map_dbl(seq(1, rv$n, by = 1), ~ input[[paste0("amt", .x)]]),
      int = purrr::map_dbl(seq(1, rv$n, by = 1), ~ input[[paste0("int", .x)]]),
      dur = purrr::map_dbl(seq(1, rv$n, by = 1), ~ input[[paste0("dur", .x)]]),
      bwt = input$bwt,
      nid = input$nid
    )
  })  # reactive
  
# Observe add button for dynamic ui
  observeEvent(input$add, {
    names <- c("amt", "int", "dur")
    index <- seq(1, rv$n, by = 1)
    purrr::walk(names, ~ { rv[[.x]][index] <- Rinput()[[.x]][index] })
    if (rv$n < 3) { rv$n <- rv$n + 1 }
  })  # observeEvent
  
# Observe remove button for dynamic ui
  observeEvent(input$rem, {
    names <- c("amt", "int", "dur")
    index <- seq(1, rv$n, by = 1)
    purrr::walk(names, ~ { rv[[.x]][index] <- Rinput()[[.x]][index] })
    if (rv$n != 1) { rv$n <- rv$n - 1 }
  })  # observeEvent
  
# Observe values of input$choose for dynamic ui
  observeEvent(input$choose, {
    rv$n <- dplyr::case_when(
      as.double(input$choose) %in% c(1, 3, 4) ~ 1,
      as.double(input$choose) == 2 ~ 2,
      as.double(input$choose) == 0 ~ rv$n)
    index <- seq(1, rv$n, by = 1)
    if (input$choose == "1") {  # "Days 1-3: 500 mg" = 1,
      rv$amt[index] <- 500
      rv$int[index] <- 24
      rv$dur[index] <- 3
    } else if (input$choose == "2") {  # "Day 1: 500mg; Days 2-5: 250 mg" = 2,
      rv$amt[index] <- c(500, 250)
      rv$int[index] <- c(24, 24)
      rv$dur[index] <- c(1, 4)
    } else if (input$choose == "3") {  # "Day 1: 1000mg" = 3,
      rv$amt[index] <- 1000
      rv$int[index] <- 24
      rv$dur[index] <- 1
    } else if (input$choose == "4") {  # "Days 1-10: 500mg" = 4
      rv$amt[index] <- 500
      rv$int[index] <- 24
      rv$dur[index] <- 10
    }  # end if
  })  # observeEvent
  
# Return dynamic ui as output
  output$regimen <- renderUI(Rui())
  
# Call module responsible for simulation and pass in inputs
# Assign output from module to rsim
  rsim <- callModule(mod_simulate_server, "sim", Rinput = Rinput, rv = rv)
  
# Return rsim object (from mod_simulate_server) to parent (app_server) 
  return(rsim)
  
}  # mod_regimen_server
