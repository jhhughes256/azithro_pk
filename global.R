# Import packages and utility functions ---------------------------------------
  library(shiny)
  library(shinydashboard)
  source("utils.R")

# Source in Shiny modules -----------------------------------------------------
  source("mod_regimen.R")
  source("mod_simulate.R")
  source("mod_plotsim.R")
  source("mod_tablesim.R")
  source("mod_infotab.R")

# Source in business logic ----------------------------------------------------
  source("fct_compile_model.R")
  source("fct_simulate_model.R")

# Set PATH for mrgsolve and run app -------------------------------------------
  if(interactive()) {
  # Get version of R
    Rversion <- paste("R-", getRversion(), sep="")
  # Set path for Rtools
    Sys.setenv(PATH = paste0(glue::glue(
      "c:/program files/R/{Rversion}/bin/x64/;", "c:/RTools/bin/;",
      "c:/RTools/mingw_64/bin/;", Sys.getenv("PATH")
    )))  # Sys.setenv
  # UI for application
    source("ui.R")
  # Server for application 
    source("server.R")
  # Run app
    runApp(shinyApp(app_ui, app_server))
  }
