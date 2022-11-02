#' Obtain GTF file of the genome annotation for a specific species from Ensembl
#'
#'
#' A function that extracts reference transcriptome from Ensembl. If a release
#' version is provided, the data will be extracted from the given version. If
#' the release version is not provided, then the default version would be the
#' latest version of Ensembl Archive.
#'
#' @param species A string indicating the specific name of the species to get
#'     reference transcriptome.
#' @param wantedVersion A double indicating the version of Ensembl archive.
#'
#'
#' @return Returns the name of the gtf.gz compressed file of the
#' reference transcriptome.
#'
#'
#' @example
#' obtainGTF("Caenorhabditis Elegans", wantedVersion=107)
#'
#' @import utils
#' @import rvest
#' @import stringr

obtainGTF <- function(species, wantedVersion=NA) {

    # species <- "CaenorhAbditis EleGans"
    species <- tolower(species)

    species <- gsub(" ", "_", species)
    # species

    if (is.na(wantedVersion)){
        # Obtain latest version
        ensemblArchives <- biomaRt::listEnsemblArchives()
        versions <- ensemblArchives$version
        wantedVersion <- suppressWarnings(na.omit(as.numeric(versions))[1])
    }

    url <- paste0("https://ftp.ensembl.org/pub/release-",
                  wantedVersion,
                  "/gtf/",
                  species)

    gtfFile <- rvest::read_html(url)
    gtfFile <- rvest::html_nodes(gtfFile, "a")
    gtfFile <- rvest::html_attr(gtfFile, "href")
    gtfFile <- stringr::str_subset(gtfFile,
                                   paste0("\\.", wantedVersion,".gtf.gz"))
    gtfFile <- gtfFile[[1]]

    # Download the file to current working directory

    utils::download.file(paste0(url,"/", gtfFile), destfile = basename(gtfFile))

    return(gtfFile)
}

