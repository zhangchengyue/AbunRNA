#' Download salmon software through bioconda.
#'
#' A function that installs Salmon software for indexing and quantification.
#' If conda is not previously installed, the function would install for use of
#' Salmon.
#'
#' @examples
#' \dontrun{
#' installSalmon()
#' }
#'
#' @references
#' Ushey K, Allaire J, Wickham H, Ritchie G (2022). rstudioapi: Safely Access
#' the RStudio API.
#' \href{https://rstudio.github.io/rstudioapi/}{Link}
#' \href{https://github.com/rstudio/rstudioapi}{Link}
#'
#' @export
#' @importFrom rstudioapi terminalSend

installSalmon <- function() {
    myTerm <- rstudioapi::terminalCreate()

    # Download conda if conda not exists
    path <- getwd()
    rstudioapi::terminalSend(myTerm, paste0("which conda > ",
                                            path, "/tmp.txt\n"))
    Sys.sleep(2)
    tmp <- read.table(file = "tmp.txt")
    tmp <- as.character(tmp)
    tmp <- paste0(tmp, collapse = " ")

    if (tmp == "conda not found") {
        # Get miniconda
        rstudioapi::terminalSend(myTerm, "brew install anaconda\n")
    }

    # Delete temporary file immediately after job done
    unlink("tmp.txt")

    # Conda config
    rstudioapi::terminalSend(myTerm,
                             "conda config --add channels conda-forge\n")
    rstudioapi::terminalSend(myTerm,
                             "conda config --add channels bioconda\n")

    # Download Salmon from conda
    rstudioapi::terminalSend(myTerm,
                             "conda create -n salmon salmon -y\n")

    # Activate salmon conda environmen
    rstudioapi::terminalSend(myTerm,
                             "conda activate salmon\n")
    return(invisible(NULL))
}

# [END]
