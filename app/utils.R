#
#' Create action button with specific colour
#' 
#' @description 
#' Creates an action button with a colour that matches shinydashboards status
#' colours
#' 
#' @param inputId The `input` slot that will be used to access the value.
#' @param label The concents of the button or link-usually a text label, but you
#' could use any other HTML, like an image.
#' @param status The shinydashboard status that will colour the button.
#' @param icon An optional `icon` to appear on the button.
#' @param width The width of the input.
#' @param style Additional style attributes to be applied to the button.
#' @param ... Named attributes to be applied to the button or link other than
#' style attributes.
#' 
#' @return An input dataset for a mrgsolve model. 
#' @keywords internal
#' @export
#' @import shiny

dashboardButton <- function(inputId, label, status, icon = NULL, width = NULL, style = "", ...) {
  if (status %in% shinydashboard:::validStatuses) {
    button_colour <- purrr::when(status,
        . == "primary" ~ "#3C8DBC",
        . == "success" ~ "#00A65A",
        . == "info" ~ "#00C0EF",
        . == "warning" ~ "#F39C12",
        . == "danger" ~ "#DD4B39"
      )
  } else {
    stop("Must specify a valid status! Please check shinydashboard documentation 
      using ?validStatuses")
  }
# Create action button to be include in ui
  actionButton(
    inputId, label, icon = icon, width = width, 
    style = paste0("border-color:#FFFFFF; background-color:", button_colour,
      "; color:#FFFFFF;", style), ...
  )  # actionButton
}

#
#' magrittr forward-pipe operator
#' 
#' @description
#' Pipe an object forward into a function or call expression. Imported from
#' magrittr because pipe is life.
#' 
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @keywords internal

`%>%` <- magrittr::`%>%`

#
#' Print `mrgsolve` model info
#' 
#' @description 
#' Prints info about the `mrgsolve` model. Works independently 
#' 
#' @param mod mrgsolve model
#' 
#' @return Model information from mrgsolve compile model object
#' @keywords internal
#' @export
#' 
print_model_info <- function(mod) {
  print(mrgsolve::param(mod))
  print("$OMEGA")
  print(mrgsolve::omat(mod))
  print("$SIGMA")
  print(mrgsolve::smat(mod))
  mrgsolve::blocks(mod, MAIN, ODE, TABLE)
}

#
#' Create list object of plot summary statistics
#' 
#' @description 
#' Generates a list of functions to be fed to a summarise_at function to provide
#' summary statistics for plotting. List consists of the median and the 
#' requested prediction intervals (percentiles).
#' 
#' @param civals Desired prediction intervals for summary statistic list
#' 
#' @return List of functions that can be used with summarise_at
#' @keywords internal
#' @export
#' 
plot_summary_fn <- function(pivals = c(0.9, 0.8, 0.6, 0.4, 0.2)) {
  purrr::map(pivals, function(pival) {
      pilo <- (1-pival)/2 # Lower percentile
      pihi <- pival+pilo # Upper percentile
      out <- list(
        as.formula(paste0("~quantile(., probs = ", pilo, ", na.rm = TRUE)")),
        as.formula(paste0("~quantile(., probs = ", pihi, ", na.rm = TRUE)"))
    )}) %>%
    unlist(recursive = FALSE) %>% 
    {purrr::list_merge(list(~median(., na.rm = TRUE)), .)} %>% 
    unlist(recursive = FALSE) %>% 
    magrittr::set_names(c("median", paste0(
      "pi", rep(pivals*100, each = 2), rep(c("lo", "hi"), times = length(pivals))
    )))  # paste0, c, set_names
}

#