#' Generate txdb object.
#'
#' A function that generates a txdb object from input quantification files,
#' which can later be used for generating count matrix of transcript abundance
#' and the downstream expression analysis pipeline.
#'
#' @param sfSeq A character vector indicating the names of salmon
#'     output .sf files.
#' @param refTrp A string indicating the path to the reference
#'     transcriptome file. The default value is NA
#' @param species A string indicating the specific name of the species to get
#'     reference transcriptome.
#' @param release A string indicating the release version of the reference.
#'     Default is the current release.
#' @param key A stirng indicating the keytype to be used.
#'
#'
#' @return Returns a txdb object
#'
#'
#' @examples
#' # Example 1:
#' # Directly obtain the latest reference transcriptome from Ensembl database
#' txi <- txdbObj(sfSeq,
#'                 species = "Caenorhabditis elegans",
#'                 key = "TXNAME")
#'
#' # Example 2:
#' # Directly obtain version 107 reference transcriptome from Ensembl database
#' txi <- txdbObj(sfSeq,
#'                 species = "Caenorhabditis elegans",
#'                 release = 107
#'                 key = "TXNAME")
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
#' Wickham H, François R, Henry L, Müller K (2022). dplyr: A Grammar of Data
#' Manipulation.
#' \href{https://dplyr.tidyverse.org}{Link}
#' \href{https://github.com/tidyverse/dplyr}{Link}
#'
#' Durinck S, Spellman P, Birney E, Huber W (2009). “Mapping identifiers for the
#'  integration of genomic datasets with the R/Bioconductor package biomaRt.”
#'  Nature Protocols, 4, 1184–1191.
#'
#' Durinck S, Moreau Y, Kasprzyk A, Davis S, De Moor B, Brazma A, Huber W
#' (2005). “BioMart and Bioconductor: a powerful link between biological
#' databases and microarray data analysis.” Bioinformatics, 21, 3439–3440.
#'
#' Lawrence M, Huber W, Pagès H, Aboyoun P, Carlson M, Gentleman R, Morgan M,
#' Carey V (2013). “Software for Computing and Annotating Genomic Ranges.” PLoS
#' Computational Biology, 9. doi: 10.1371/journal.pcbi.1003118,
#' \href{http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003118}{Link}.
#'
#'
#'
#' @export
#' @import utils
#' @import rvest
#' @import stringr
#' @import biomaRt
#' @import dplyr
#' @import GenomicFeatures
#' @import tximport
#' @import AnnotationDbi


txdbObj <- function(sfSeq,
                    refTrp = NA,
                    species=NA,
                    release=NA,
                    key = "TXNAME"){

    if (typeof(sfSeq) != "character") {
        stop("Please provide valid .sf quantification files.")
    }
    # Make a txdb object
    # If no reference transcriptome provided, obtain the annotation from Ensembl
    if (is.na(refTrp)) {

        # If no species name is provided, stop with an error message
        if (is.na(species)) {

            stop("No valid reference transcriptome file or species provided.")
        }

        refTrp <- obtainGTF(species, release)

        # Unzip the gtf file
        system(paste("gunzip", refTrp))

        # Reformat the file name used for further process
        refTrp <- sub(".gz", "", refTrp)

    }
    txdb <- GenomicFeatures::makeTxDbFromGFF(refTrp)

    # Obtain txname keys
    k <- biomaRt::keys(txdb, keytypes = key)

    txGene <- AnnotationDbi::select(txdb, keys = k, columns = "TXNAME", keytype = "GENEID")

    txGene <- txGene[, c("TXNAME", "GENEID")]

    txi <- tximport::tximport(sfSeq, type="salmon", tx2gene=txGene)

    return(txi)
}

# [END]
