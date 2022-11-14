#' Download salmon software through bioconda.
#'
#'
#' A function that downloads salmon for indexing and quantification if required.
#' If conda is not previously installed, the function would install for use of
#' salmon.
#'
#'
#' @example
#' install
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

install <- function(){
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

