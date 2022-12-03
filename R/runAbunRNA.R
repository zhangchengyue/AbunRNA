#' Launch Shiny App for AbunRNA
#'
#' A function that launches the Shiny app for AbunRNA。
#' The shiny app permits to generate count matrix, plot heatmap and
#' generate PCA plot.
#'
#'
#' @return No return value but open up a Shiny page.
#'
#' @examples
#' \dontrun{
#'
#' AbunRNA::runAbunRNA()
#' }
#'
#' @references
#' Grolemund, G. (2015). Learn Shiny - Video Tutorials. \href{https://shiny.rstudio.com/tutorial/}{Link}
#'
#' @export
#' @importFrom shiny runApp

runAbunRNA <- function() {
    appDir <- system.file("shiny-scripts",
                          package = "AbunRNA")
    shiny::runApp(appDir, display.mode = "normal")
    return()
}
# [END]
