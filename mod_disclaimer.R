# Module UI
  
#' @title   mod_disclaimer_ui and mod_disclaimer_server
#' @description  Shiny module which handles the modal dialogue box which
#'   displays the disclaimer requiring acknowledgement prior to use of the
#'   application by the user.
#'
#' @param id shiny id
#' @param input internal
#' @param output internal
#' @param session internal
#'
#' @rdname mod_disclaimer
#'
#' @keywords internal
#' @export 
#' @import shiny

mod_disclaimer_ui <- function(id) {
# Define namespace function for IDs
  ns <- NS(id)
# Create tagList to be used in the UI
  modalDialog(
    title = fluidRow(
      div(
        icon("exclamation-circle"), 
        style = "font-size:80px; color:#F39C12"
      ),
      div(h1("ATTENTION"), align = "center"),
      align = "center"
    ),
    tagList(
      fluidRow(
        column(12,
          p(paste("This application was developed to aid dosing regimen selection", 
            "for azithromycin in the treatment of COVID-19 conducted under the", 
            "controlled environment of a registered clinical trial.")),
          p(paste("Azithromycin is not approved by the FDA for treatment of", 
            "COVID-19. No recommendations are made for the treatment of COVID-19",
            "with azithromycin outside of a registered clinical trial.")),
          p(paste("While this web application allows the comparison of in-vitro", 
            "effective concentration values with tissue exposure, these values may" , 
            "not be indicative of effective concentrations in patients with COVID-19.")),
          p(paste("The published model used for the simulation of azithromycin",
            "treatment regimens is documented in the manuscript below."), 
            tags$blockquote(paste("Hughes, J.H., Sweeney, K., Ahadieh, S. and",
              "Ouellet, D. (2020), Predictions of Systemic, Intracellular, and", 
              "Lung Concentrations of Azithromycin with Different Dosing Regimens",
              "used in COVID-19 Clinical Trials. CPT Pharmacometrics Syst. Pharmacol."), 
              a("doi:10.1002/psp4.12537", 
                href = "https://ascpt.onlinelibrary.wiley.com/doi/abs/10.1002/psp4.12537"),
              style = "font-size:14px"
          ))
        ),
        align = "justified"
      )
    ),
    footer = div(modalButton("I acknowledge the information above"), align = "center"),
    size = "l", easyClose = FALSE
  )
}  # mod_disclaimer_ui