# Module UI
  
#' @title   mod_tablesim_ui and mod_tablesim_server
#' @description  Shiny module which handles display of dosing regimen text and
#' table of metrics. Includes data processing required to generate the table and 
#' regimen strings. Child module of mod_plotsim.
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
      uiOutput(ns("regout")),
      tableOutput(ns("tabout"))
    ),  # div
    div(class = "col-sm-12 col-md-6", align = "center",
      uiOutput(ns("regsave")),
      tableOutput(ns("tabsave"))
    ),  # div
    column(12,
      p(paste("Abbreviations: AM - alveolar macrophage; PI - prediction interval;", 
      "EC50: 50% effective concentration; EC90: 90% effective concentration;",
      "WBC: white blood cell predictions for polymorphonuclear leukocytes,",
      "mononuclear leukocytes and peripheral blood monocytes."))
    )  # column
  )  # tagList
}  # mod_tablesim_ui

# Module Server
    
#' @rdname mod_tablesim
#' @export
#' @keywords internal

mod_tablesim_server <- function(input, output, session, rsim) {
# Define namespace function for IDs
  ns <- session$ns
  
# Define function for generating dose regimen text
# * If there is a regimen available for a given output
#     + Create a list of strings that describe the dosing regimen based on input
#     + Convert to HTML object with strong and h5 tags
#     + When more sets of HTML would be created by the other output
#         - Add additional br tags to ensure tables are level
# * Otherwise if there is a regimen available for the other output
#     + return br tags equal to the sets of HTML that would be created by the 
#       other output
# * Otherwise there are no regimens available, return a blank string
  fct_reggen <- function(reg, other) {
    if (!is.null(reg)) {
      seq_along(reg$amt) %>%
        purrr::map(function(i) {
          dur_txt <- paste0("Day", ifelse(reg$dur[[i]] > 1, "s ", " "))
          dur_val <- reg$dur %>% 
            purrr::when(
              .[[i]] == 1 & i == 1 ~ "1",
              .[[i]] != 1 & i == 1 ~ paste("1", .[[i]], sep = "-"),
              .[[i]] == 1 & i != 1 ~ paste(.[[i-1]] + 1),
              .[[i]] != 1 & i != 1 ~ paste(.[[i-1]] + 1, .[[i-1]] + .[[i]], sep = "-")
            )  # when
          dur_str <- paste0(dur_txt, dur_val, ":")
          amt_str <- paste(reg$amt[[i]], "mg")
          int_str <- paste("every", reg$int[[i]], "hours")
          paste(dur_str, amt_str, int_str)
        }) %>%  # map
        purrr::map(strong) %>%
        purrr::map(h5) %>%
        purrr::when(
          is.null(other) ~ .,
          length(.) >= length(other$amt) ~ .,
          length(.) < length(other$amt) ~ list(.,
            purrr::map(
              seq(length(reg$amt) + 1, length(other$amt), by = 1), ~ h5(br()))) %>%
            unlist(recursive = FALSE)
        ) %>%  # when
        return()
    } else if (!is.null(other)) {
      seq_along(other$amt) %>%
        purrr::map(~ h5(br())) %>%
        return()
    } else {
      return("")
    }
  }
  
# Render dose regimen text to be displayed in the application UI
  output$regout <- renderUI({
    fct_reggen(rsim$out_reg, rsim$save_reg)
  })
  output$regsave <- renderUI({
    fct_reggen(rsim$save_reg, rsim$out_reg)
  })
  
# Define function for generating application tables
# * Determine metric label based on number of indiviudals simulated
# * When mrgsolve output is available
#     + Filter for the final time record of each individual
#     + Convert time above EC to days and round to one decimal
#     + Determine the median of all columns (for simulation of multiple IDs)
#     + Produce a table for display
# * Otherwise when mrgsolve output is NULL produce an empty table
  fct_tabgen <- function(rsim, label) {
    nid <- length(unique(rsim$ID))
    metric <- ifelse(nid > 1, "Median (90% PI)", "Value")
    rsim %>%
      purrr::when(
        !is.null(.) ~ rsim %>%
          dplyr::filter(time == max(time)) %>%
          dplyr::select(tidyselect::contains("EC")) %>%
          dplyr::mutate_at(dplyr::vars(tidyselect::contains("TEC")), 
            magrittr::divide_by, 24) %>%
          dplyr::group_by(AZEC50, AZEC90) %>%
          dplyr::summarise_all(plot_summary_fn(ifelse(nid > 1, 0.9, NA))) %>%
          purrr::when(
            nid != 1 ~ tidyr::pivot_longer(., tidyselect::contains("TEC")) %>%
              tidyr::separate("name", c("metric", "stat"), sep = "_") %>%
              dplyr::mutate(value = round(value, 1)) %>%
              tidyr::pivot_wider(names_from = "stat", values_from = "value") %>%
              dplyr::mutate(value = paste0(median, " (", pi90lo, " - ", pi90hi, ")")) %>%
              dplyr::select(-median, -pi90lo, -pi90hi) %>%
              tidyr::pivot_wider(names_from = "metric", values_from = "value"),
            nid == 1 ~ mutate_at(., dplyr::vars(-dplyr::group_cols()), round, 1)) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(NID = nid) %>%
          dplyr::select(
            "Number of Individuals" = NID,
            "EC50 (ng/mL)" = AZEC50,
            "Time above EC50 in AM (days)" = ALMATEC50,
            "Time above EC50 in lung (days)" = LUNGTEC50,
            "EC90 (ng/mL)" = AZEC90,
            "Time above EC90 in AM (days)" = ALMATEC90,
            "Time above EC90 in lung (days)" = LUNGTEC90) %>%
          dplyr::mutate_all(as.character),
        is.null(.) ~ tibble::tibble(
          "Number of Individuals" = "-",
          "EC50 (ng/mL)" = "-",
          "Time above EC50 in AM (days)" = "-",
          "Time above EC50 in lung (days)" = "-",
          "EC90 (ng/mL)" = "-",
          "Time above EC90 in AM (days)" = "-",
          "Time above EC90 in lung (days)" = "-")) %>%
      tidyr::pivot_longer(cols = tidyselect::everything(), 
        names_to = label, values_to = metric)
  }

# Render tables to be displayed in the application
  output$tabout <- renderTable({
    fct_tabgen(rsim$out, "Current Regimen")
  })
  output$tabsave <- renderTable({
    fct_tabgen(rsim$save, "Saved Regimen")
  })

}  # mod_tablesim_server
