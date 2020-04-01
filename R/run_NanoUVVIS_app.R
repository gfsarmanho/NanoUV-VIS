#' Shiny UI
#' 
#' \code{NanoUVVIS.App} runs the  Shiny web application.
#' 
#' @examples
#' 
#' #NanoUVVIS.App()
#'

#' @export
run_NanoUVVIS_app <- function() {
  appDir = system.file("App", package = "NanoUVVIS.App")
  if (appDir == "") {
    stop("Could not find app. Try re-installing `NanoUVVIS.App`.", call. = FALSE)
  }
  shiny::runApp(appDir, launch.browser = TRUE, display.mode = "normal")
}