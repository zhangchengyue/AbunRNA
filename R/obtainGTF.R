#' Obtain GTF file of the genome annotation for a specific species from Ensembl
#'
#' A function that extracts reference transcriptome from Ensembl. If a release
#' version is provided, the data will be extracted from the given version. If
#' the release version is not provided, then the default version would be the
#' latest version of Ensembl Archive.
#'
#' @param species A string indicating the specific name of the species to get
#'     reference transcriptome.
#' @param wantedVersion A double indicating the version of Ensembl archive.
#' @param download A boolean indicating whether to download the file or not.
#'
#' @return Returns the name of the gtf.gz compressed file of the
#' reference transcriptome.
#'
#' @example
#' \dontrun{
#' obtainGTF("Caenorhabditis Elegans", wantedVersion=107)
#' }
#'
#' @references
#' Ensembl 2022, Nucleic Acids Research, Volume 50, Issue D1,
#' 7 January 2022, Pages D988–D995,
#' \href{https://doi.org/10.1093/nar/gkab1049}{Link}
#'
#' Wickham H (2022). rvest: Easily Harvest (Scrape) Web Pages.
#' \href{https://rvest.tidyverse.org/}
#' \href{https://github.com/tidyverse/rvest}
#'
#' Wickham H (2022). stringr: Simple, Consistent Wrappers for
#' Common String Operations.
#' \href{http://stringr.tidyverse.org}
#' \href{https://github.com/tidyverse/stringr}
#'
#' @import utils
#' @import rvest
#' @import stringr

obtainGTF <- function(species = NA, wantedVersion = NA, download = F) {
    if (is.na(species)) {
        stop("Species name not provided.")
    } else {
        ;
    }

    # Format the input name of species
    species <- tolower(species)

    species <- gsub(" ", "_", species)

    if (is.na(wantedVersion)) {
        # Obtain latest version
        ensemblArchives <- biomaRt::listEnsemblArchives()
        versions <- ensemblArchives$version
        versions <- as.numeric(versions)
        versions <- na.omit(versions)
        wantedVersion <- suppressWarnings(versions[1])
    } else {
        ;
    }

    url <- paste0("https://ftp.ensembl.org/pub/release-",
                  wantedVersion,
                  "/gtf/",
                  species)

    code <- tryCatch({gtfFile <- rvest::read_html(url)},
                     error = function(e) {return(1)})
    if (typeof(code) == "double") {
        stop("Species name cannot be recognized. Please try again.")
    } else {
        ;
    }
    gtfFile <- rvest::html_nodes(gtfFile, "a")
    gtfFile <- rvest::html_attr(gtfFile, "href")
    gtfFile <- stringr::str_subset(gtfFile,
                                   paste0("\\.", wantedVersion,".gtf.gz"))
    gtfFile <- gtfFile[[1]]

    # Download the file to current working directory
    if (download == T) {
        utils::download.file(paste0(url,"/", gtfFile),
                             destfile = basename(gtfFile))
    } else {
        ;
    }

    return(gtfFile)
}

# [END]
