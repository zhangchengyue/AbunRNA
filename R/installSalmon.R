#' Download salmon software through bioconda.
#'
#'
#' A function that installs Salmon software for indexing and quantification.
#' If conda is not previously installed, the function would install for use of
#' Salmon.
#'
#'
#' @example
#' install()
#'
#'
#' @references
#' Ushey K, Allaire J, Wickham H, Ritchie G (2022). rstudioapi: Safely Access
#' the RStudio API.
#' \href{https://rstudio.github.io/rstudioapi/}
#' \href{https://github.com/rstudio/rstudioapi}
#'
#'
#'
#' @export
#' @import rstudioapi

installSalmon <- function(){
    myTerm <- rstudioapi::terminalCreate()

    # Download conda if conda not exists
    rstudioapi::terminalSend(myTerm, "which conda > tmp.txt\n")
    tmp <- paste0(as.character(read.table(file = "tmp.txt")), collapse = " ")
    if (tmp == "conda not found"){
        # get miniconda
        rstudioapi::terminalSend(myTerm, "brew install anaconda\n")
    }
    unlink("tmp.txt")
    rstudioapi::terminalSend(myTerm,
                             "conda config --add channels conda-forge\n")
    rstudioapi::terminalSend(myTerm,
                             "conda config --add channels bioconda\n")

    # Download salmon from conda
    rstudioapi::terminalSend(myTerm,
                             "conda create -n salmon salmon -y\n")

    # Activate salmon conda environment
    rstudioapi::terminalSend(myTerm,
                             "conda activate salmon\n")
}

# [END]
