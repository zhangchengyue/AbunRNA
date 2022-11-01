#' Generate a count matrix from salmon output files.
#'
#' A function that process the input salmon output .sf files, extract genes from
#' the provided gene annotation of the species, and generate a count
#' matrix to record the expression abundance of the genes for raw sequences.
#'
#'
#' @param sfSeq A character vector indicating the names of salmon
#'     output .sf files.
#' @param refTrp A string indicating the path to the reference
#'     transcriptome file. The default value is NA
#' @param sampleNames A character vector indicating sample names of the raw
#'     sequences.
#' @param type A string indicating the tool used for generating transcript
#'     expression abundance from raw sequence. The default tool is set to salmon.
#' @param species A string indicating the specific name of the species to get
#'     reference transcriptome.
#' @param release A string indicating the release version of the reference.
#'     Default is the current release.
#' @param abunCSV A character string for the name of the output abundance CSV
#'
#'
#' @return Returns a data frame where each row indicates a gene, each column
#' indicates a sample, each cell indicates the expression abundance of a gene
#' in a sample.
#'
#'
#' @example
#' # Provide the specific name of the species to directly obtain the
#' # reference transcriptome from Ensembl database. Database would be specific
#' # to a certian version if release version is provided (e.g. 107 in the
#' # following example).
#'
#' # Access .sf files generated from salmon that are available in this package
#' cw1_quants <- system.file("extdata", "cw1_quants", "quant.sf", package = "AbunRNA")
#' cl1_quants <- system.file("extdata", "cl1_quants", "quant.sf", package = "AbunRNA")
#' cg1_quants <- system.file("extdata", "cg1_quants", "quant.sf", package = "AbunRNA")
#' sfSe <- c(cw1_quants, cl1_quants, cg1_quants)
#'
#' # Name the samples correspondingly
#' samples <- c("WT_WC_1", "lf_WC_1", "gf_WC_1")
#' output <- generateMatrix(sfSeq = sfSe,
#'                       species = "Caenorhabditis elegans",
#'                       release = 107,
#'                       sampleNames = samples,
#'                       abunCSV = "abunOUT")
#' View(output)
#'
#'
#'
#' @references
#' Ensembl 2022, Nucleic Acids Research, Volume 50, Issue D1,
#' 7 January 2022, Pages D988–D995,
#' \href{https://doi.org/10.1093/nar/gkab1049}{Link}
#'
#' Durinck S, Spellman P, Birney E, Huber W (2009). “Mapping identifiers for the
#'  integration of genomic datasets with the R/Bioconductor package biomaRt.”
#'  Nature Protocols, 4, 1184–1191.

#' Durinck S, Moreau Y, Kasprzyk A, Davis S, De Moor B, Brazma A, Huber W
#' (2005). “BioMart and Bioconductor: a powerful link between biological
#' databases and microarray data analysis.” Bioinformatics, 21, 3439–3440.
#'
#' Wickham H, François R, Henry L, Müller K (2022). dplyr: A Grammar of Data
#' Manipulation.
#' \href{https://dplyr.tidyverse.org}{Link}
#' \href{https://github.com/tidyverse/dplyr}{Link}
#'
#' Müller K, Wickham H (2022). tibble: Simple Data Frames.
#' \href{https://tibble.tidyverse.org/}{Link}
#' \href{https://github.com/tidyverse/tibble}{Link}
#'
#' Lawrence M, Huber W, Pagès H, Aboyoun P, Carlson M, Gentleman R, Morgan M,
#' Carey V (2013). “Software for Computing and Annotating Genomic Ranges.” PLoS
#' Computational Biology, 9. doi: 10.1371/journal.pcbi.1003118,
#' \href{http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003118}{Link}.
#'
#'
#'
#' @export
#' @import biomaRt
#' @import dplyr
#' @import tibble
#' @import GenomicFeatures
#' @import tximport

library(GenomicFeatures)
generateMatrix <- function (sfSeq,
                            refTrp=NA,
                            sampleNames,
                            type="salmon",
                            species=NA,
                            release=NA,
                            abunCSV="abunCSV") {

    if (typeof(sfSeq) != "character") {
        stop("Please provide valid .sf quantification files.")
    }


    # Make a txdb object for the provided species.
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

    # Make a txdb object
    txdb <- GenomicFeatures::makeTxDbFromGFF(refTrp)

    # Obtain txname keys
    k <- biomaRt::keys(txdb, keytypes = "TXNAME")

    txGene <- select(txdb, keys = k, columns = "TXNAME", keytype = "GENEID")

    txGene <- txGene[, c("TXNAME", "GENEID")]

    txi <- tximport::tximport(sfSeq, type="salmon", tx2gene=txGene)

    abunFile <- paste(abunCSV, ".csv")

    write.csv(txi$abundance, file = abunFile)

    abundance <- read.csv(abunFile)

    abunAnnot <- data.frame(abundance)

    colnames(abunAnnot) <- c("Gene ID", sampleNames)

    # A function used as a parameter of the "apply" function to filter out the
    # unexpressed genes.
    unEx <- function(x) {
        return(!all(x==0))
    }

    # Get rid of unexpressed genes
    abunAnnot[apply(abunAnnot, 1, unEx), ]

    return(abunAnnot)
}


#' Obtain GTF file of a specific species from Ensmbl latest version
#'
#'
#' Helper function for generateMatrix to obtain reference transcriptome from
#' the latest version of Ensembl.
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

    gtfFile <- rvest::read_html(url) %>%
        html_nodes("a") %>%
        html_attr("href") %>%
        stringr::str_subset(paste0("\\.",wantedVersion,".gtf.gz")) %>%
        .[[1]]

    # Download the file to current working directory

    utils::download.file(paste0(url,"/", gtfFile), destfile = basename(gtfFile))
    return(gtfFile)
}
