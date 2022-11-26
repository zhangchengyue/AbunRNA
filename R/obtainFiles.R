#' Obtain complete cDNA fa.gz file for a specific species from Ensembl
#'
#' A function that extracts complete cDNA fa.gz file from Ensembl. If a release
#' version is provided, the data will be extracted from the given version. If
#' the release version is not provided, then the default version would be the
#' latest version of Ensembl Archive.
#'
#' @param species A string indicating the specific name of the species.
#'     The default value is NA.
#' @param wantedVersion A double indicating the version of Ensembl archive.
#'     The default value is NA.
#' @param download A boolean indicating whether to download the file or not.
#'     The default value is FALSE.
#'
#'
#' @return Returns the name of the fa.gz compressed file of cDNA.
#'
#'
#' @examples
#' \dontrun{
#' obtainCDNA(species = "Caenorhabditis Elegans", wantedVersion = 107)
#' }
#'
#'
#' @references
#' Ensembl 2022, Nucleic Acids Research, Volume 50, Issue D1,
#' 7 January 2022, Pages D988–D995,
#' \href{https://doi.org/10.1093/nar/gkab1049}{Link}
#'
#' Wickham H (2022). rvest: Easily Harvest (Scrape) Web Pages.
#' \href{ https://CRAN.R-project.org/package=rvest}{Link}
#'
#' Wickham H (2022). stringr: Simple, Consistent Wrappers for
#' Common String Operations.
#' \href{https://CRAN.R-project.org/package=stringr}{Link}
#'
#' @export
#' @importFrom rvest read_html
#' @importFrom rvest html_nodes
#' @importFrom rvest html_attr
#' @importFrom stringr str_subset
#' @importFrom utils download.file
#' @importFrom biomaRt listEnsemblArchives


obtainCDNA <- function(species = NA, wantedVersion = NA, download = F) {

    if (is.na(species)) {
        stop("Species name not provided.")
    } else {
        ;
    }

    # Format input species name
    species <- tolower(species)

    species <- gsub(" ", "_", species)

    if (is.na(wantedVersion)) {
        # Obtain latest version
        ensemblArchives <- biomaRt::listEnsemblArchives()
        versions <- ensemblArchives$version
        versions <- suppressWarnings(as.numeric(versions))
        versions <- suppressWarnings(na.omit(versions))
        wantedVersion <- suppressWarnings(versions[1])
    } else {
        ; # Wanted version provided. Do nothing
    }
    url <- paste0("https://ftp.ensembl.org/pub/release-",
                  wantedVersion,
                  "/fasta/",
                  species,
                  "/cdna/")
    code <- tryCatch({cdna <- rvest::read_html(url)},
             error = function(e) {return(1)})
    if (typeof(code) == "double") {
        stop("Species name cannot be recognized. Please try again.")
    } else {
        ;
    }
    cdna <- rvest::html_nodes(cdna, "a")
    cdna <- rvest::html_attr(cdna, "href")
    cdna <- stringr::str_subset(cdna, ".*\\.cdna.all.fa.gz$")
    cdna <- cdna[[1]]

    # Download the file to current working directory
    if (download == T) {
        utils::download.file(paste0(url,"/", cdna),
                             destfile = basename(cdna))
    } else {
        ;
    }

    return(cdna)
}





#' Obtain DNA toplevel fa.gz file for a specific species from Ensembl
#'
#' A function that extracts DNA toplevel fa.gz file from Ensembl. If a release
#' version is provided, the data will be extracted from the given version. If
#' the release version is not provided, then the default version would be the
#' latest version of Ensembl Archive.
#'
#' @param species A string indicating the specific name of the species.
#'     The default value is NA.
#' @param wantedVersion A double indicating the version of Ensembl archive.
#'     The default value is NA.
#' @param download A boolean indicating whether to download the file or not.
#'     The default value is FALSE.
#'
#' @return Returns the name of the fa.gz compressed file of DNA.
#'
#' @examples
#' \dontrun{
#' obtainDNA("Caenorhabditis Elegans", wantedVersion=107)
#' }
#'
#' @references
#' Ensembl 2022, Nucleic Acids Research, Volume 50, Issue D1,
#' 7 January 2022, Pages D988–D995,
#' \href{https://doi.org/10.1093/nar/gkab1049}{Link}
#'
#' Wickham H (2022). rvest: Easily Harvest (Scrape) Web Pages.
#' \href{https://rvest.tidyverse.org/}{Link}
#' \href{https://github.com/tidyverse/rvest}{Link}
#'
#' Wickham H (2022). stringr: Simple, Consistent Wrappers for
#' Common String Operations.
#' \href{http://stringr.tidyverse.org}{Link}
#' \href{https://github.com/tidyverse/stringr}{Link}
#'
#' @export
#' @import utils
#' @importFrom rvest read_html
#' @importFrom rvest html_nodes
#' @importFrom rvest html_attr
#' @importFrom stringr str_subset
#' @importFrom biomaRt listEnsemblArchives


obtainDNA <- function(species = NA, wantedVersion = NA, download = F) {
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
        versions <- suppressWarnings(as.numeric(versions))
        versions <- suppressWarnings(na.omit(versions))
        wantedVersion <- suppressWarnings(versions[1])
    } else {
        ;
    }
    url <- paste0("https://ftp.ensembl.org/pub/release-",
                  wantedVersion,
                  "/fasta/",
                  species,
                  "/dna/")

    code <- tryCatch({dna <- rvest::read_html(url)},
                     error = function(e) {return(1)})
    if (typeof(code) == "double") {
        stop("Species name cannot be recognized. Please try again.")
    } else {
        ;
    }
    dna <- rvest::html_nodes(dna, "a")
    dna <- rvest::html_attr(dna, "href")
    dna <- stringr::str_subset(dna, ".*\\.dna.toplevel.fa.gz$")
    dna <- dna[[1]]

    # Download the file to current working directory
    if (download == T) {
        utils::download.file(paste0(url, "/", dna),
                             destfile = basename(dna))
    } else {
        ;
    }

    return(dna)
}






#' Obtain GTF file of the genome annotation for a specific species from Ensembl
#'
#' A function that extracts reference transcriptome from Ensembl. If a release
#' version is provided, the data will be extracted from the given version. If
#' the release version is not provided, then the default version would be the
#' latest version of Ensembl Archive.
#'
#' @param species A string indicating the specific name of the species to get
#'     reference transcriptome. The default value is NA.
#' @param wantedVersion A double indicating the version of Ensembl archive.
#'     The default value is NA.
#' @param download A boolean indicating whether to download the file or not.
#'     The default value is FALSE.
#'
#' @return Returns the name of the gtf.gz compressed file of the
#' reference transcriptome.
#'
#' @examples
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
#' \href{https://rvest.tidyverse.org/}{Link}
#' \href{https://github.com/tidyverse/rvest}{Link}
#'
#' Wickham H (2022). stringr: Simple, Consistent Wrappers for
#' Common String Operations.
#' \href{http://stringr.tidyverse.org}{Link}
#' \href{https://github.com/tidyverse/stringr}{Link}
#'
#' @export
#' @importFrom rvest read_html
#' @importFrom rvest html_nodes
#' @importFrom rvest html_attr
#' @importFrom stringr str_subset
#' @importFrom utils download.file
#' @importFrom biomaRt listEnsemblArchives

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
        versions <- suppressWarnings(as.numeric(versions))
        versions <- suppressWarnings(na.omit(versions))
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
