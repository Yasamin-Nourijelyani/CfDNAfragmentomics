# Purpose: Allow for use of shiny
# Author: Yasamin Nouri Jelyani
# Date: 2022-12-05
# Version: 0.1.0
# Bugs and Issues:

#------------------------Exported Function------------------------



#'Use the shiny CfDNAfragmentomics
#'
#' Launch Shiny App for CfDNAfragmentomics
#'
#' A function that launches the Shiny app for CfDNAfragmentomics
#' The purpose of this app is only to run the Shiny
#' app. The code has been placed in \code{./inst/shiny-scripts}.
#'
#' @return No return value but open up a Shiny page.
#'
#' @examples
#' \dontrun{
#'}
#'
#' @references
#' 1- Grolemund, G. (2015). Learn Shiny - Video Tutorials.
#'  \href{https://shiny.rstudio.com/tutorial/}{Link}
#'
#' 2- Silva, A. (2022) TestingPackage: An Example R Package For BCB410H.
#' Unpublished. URL https://github.com/anjalisilva/TestingPackage.
#'
#' @export
#' @importFrom shiny runApp
#' @importFrom shinyalert shinyalert

runShinyCfDNAfragmentomics <- function() {
  shinyalert::useShinyalert
  #read shiny scripts
  appDir <- system.file("shiny-scripts",
                        package = "CfDNAfragmentomics")
  actionShiny <- shiny::runApp(appDir, display.mode = "normal")
  return(actionShiny)
}
# [END]
