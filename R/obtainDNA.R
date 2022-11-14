#' Obtain DNA toplevel fa.gz file for a specific species from Ensembl
#'
#'
#' A function that extracts DNA toplevel fa.gz file from Ensembl. If a release
#' version is provided, the data will be extracted from the given version. If
#' the release version is not provided, then the default version would be the
#' latest version of Ensembl Archive.
#'
#' @param species A string indicating the specific name of the species.
#' @param wantedVersion A double indicating the version of Ensembl archive.
#' @param download A boolean indicating whether to download the file or not.
#'
#'
#' @return Returns the name of the fa.gz compressed file of DNA.
#'
#'
#' @example
#' obtainDNA("Caenorhabditis Elegans", wantedVersion=107)
#'
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
#'
#' @import utils
#' @import rvest
#' @import stringr


obtainDNA <- function(species, wantedVersion=NA, download = F) {

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
    wantedVersion = 105
    url <- paste0("https://ftp.ensembl.org/pub/release-",
                  wantedVersion,
                  "/fasta/",
                  species,
                  "/dna/")
    url
    dna <- rvest::read_html(url)
    dna <- rvest::html_nodes(dna, "a")
    dna <- rvest::html_attr(dna, "href")
    dna <- stringr::str_subset(dna, ".*\\.dna.toplevel.fa.gz$")
    dna <- dna[[1]]

    # Download the file to current working directory
    if (download == T){
        utils::download.file(paste0(url,"/", dna),
                             destfile = basename(dna))}

    return(dna)
}

# [END]
