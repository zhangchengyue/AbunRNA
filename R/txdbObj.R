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
#' @return Returns a txdb object
#'
#' @examples
#' cw1Quants <- system.file("extdata",
#'                           "cw1_quants",
#'                           "quant.sf",
#'                           package = "AbunRNA")
#' cl1Quants <- system.file("extdata",
#'                           "cl1_quants",
#'                           "quant.sf",
#'                           package = "AbunRNA")
#' cg1Quants <- system.file("extdata",
#'                           "cg1_quants",
#'                           "quant.sf",
#'                           package = "AbunRNA")
#' sfSe <- c(cw1Quants, cl1Quants, cg1Quants)
#' # Example 1:
#' # Directly obtain the latest reference transcriptome from Ensembl database
#' \dontrun{
#' txi <- txdbObj(sfSeq = sfSe,
#'                 species = "Caenorhabditis elegans",
#'                 key = "TXNAME")
#' }
#'
#' # Example 2:
#' # Directly obtain version 107 reference transcriptome from Ensembl database
#' \dontrun{
#' txi <- txdbObj(sfSeq = sfSe,
#'                 species = "Caenorhabditis elegans",
#'                 release = 107,
#'                 key = "TXNAME")
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
#'
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
#'
#' @export
#' @import utils
#' @importFrom tximport tximport
#' @importFrom biomaRt keys
#' @importFrom biomaRt keytypes
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom AnnotationDbi select
#' @importFrom stringr str_subset
#' @importFrom biomaRt listEnsemblArchives
#' @importFrom rvest read_html
#' @importFrom rvest html_nodes
#' @importFrom rvest html_attr


txdbObj <- function(sfSeq = NA,
                    refTrp = NA,
                    species = NA,
                    release = NA,
                    key = "TXNAME") {

    if (typeof(sfSeq) != "character") {
        stop("Please provide valid .sf quantification files.")
    } else {
        ;
    }

    # Make a txdb object
    # If no reference transcriptome provided, obtain the annotation from Ensembl
    if (is.na(refTrp)) {

        # If no species name is provided, stop with an error message
        if (is.na(species)) {

            stop("No valid reference transcriptome file or species provided.")
        } else {
            ;
        }

        refTrp <- obtainGTF(species, wantedVersion = release, download = T)

        # Unzip the gtf file
        system(paste("gunzip", refTrp))

        # Reformat the file name used for further process
        refTrp <- sub(".gz", "", refTrp)

    } else {
        ;
    }
    txdb <- GenomicFeatures::makeTxDbFromGFF(refTrp)

    # Obtain txname keys
    allKeys <- biomaRt::keytypes(txdb)
    if (! (key %in% allKeys)) {
        stop(paste("Key doesn't exist. Available keytypes: ",
                   allKeys, sep = ",", collapse = ", "))
    } else {
        ;
    }
    k <- biomaRt::keys(txdb, keytypes = key)

    txGene <- AnnotationDbi::select(txdb, keys = k,
                                    columns = "TXNAME",
                                    keytype = "GENEID")

    txGene <- txGene[ , c("TXNAME", "GENEID")]

    txi <- tximport::tximport(sfSeq, type="salmon", tx2gene=txGene)

    return(txi)
}

# [END]
