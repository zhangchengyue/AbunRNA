#' Generate a count matrix and Plot heatmap
#'
#' A function that process the input salmon output .sf files, extract genes from
#' the provided gene annotation of the species, and generate a count
#' matrix that records the expression abundance of the genes. If indicated, the
#' function can also plot a heatmap to visualize the transcript abundance among
#' wild type and mutant samples.
#'
#' @param sfSeq A character vector indicating the names of salmon
#'     output .sf files. The default value is NA.
#' @param refTrp A string indicating the path to the reference
#'     transcriptome file. The default value is NA.
#' @param sampleNames A character vector indicating sample names of the raw
#'     sequences. The default value is NA.
#' @param type A string indicating the tool used for generating transcript
#'     expression abundance from raw sequence. The default tool is salmon.
#' @param keytype A stirng indicating the keytype to be used. Default
#'     setting is TXNAME.
#' @param element A string indicating which element should be generated from the
#'     matrix. Elements include: "abundance", "count", "length". The default
#'     value is set to "abundance".
#' @param species A string indicating the specific name of the species to get
#'     reference transcriptome. The default value is NA.
#' @param release A string indicating the release version of the reference.
#'     Default is the current release. The default value is NA.
#' @param outputCSV A boolean indicating whether to output the matrix as CSV.
#'     The default value is FALSE.
#' @param nameCSV A character string indicating the basename of the output CSV.
#'     The default output file name is "matrixCSV".
#' @param heatmap A boolean indicating whether to plot the heat map. The default
#'     value is TRUE.
#' @param head A boolean indicating whether to use only the first 6 lines of
#'     the count matrix. This is useful to view as an example when the matrix is
#'     large. The default value is TRUE.
#'
#' @return Returns a data frame where each row indicates a gene, each column
#' indicates a sample, each cell indicates the expression abundance of a gene
#' in a sample. This can be used for downstream analysis.
#'
#' @examples
#' # Provide the specific name of the species to directly obtain the
#' # reference transcriptome from Ensembl database. Database would be specific
#' # to a certian version if release version is provided (e.g. 107 in the
#' # following example).
#'
#' # Access .sf files generated from salmon that are available in this package
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
#'
#' # Name the samples correspondingly
#' samples <- c("WT_WC_1", "lf_WC_1", "gf_WC_1")
#' output <- generateMatrix(sfSeq = sfSe,
#'                       species = "Caenorhabditis elegans",
#'                       release = 107,
#'                       sampleNames = samples,
#'                       outputCSV = FALSE)
#'
#' @references
#' Ensembl 2022, Nucleic Acids Research, Volume 50, Issue D1,
#' 7 January 2022, Pages D988–D995,
#' \href{https://doi.org/10.1093/nar/gkab1049}{Link}
#'
#' Durinck S, Spellman P, Birney E, Huber W (2009). “Mapping identifiers for the
#'  integration of genomic datasets with the R/Bioconductor package biomaRt.”
#'  Nature Protocols, 4, 1184–1191.
#'
#' Durinck S, Moreau Y, Kasprzyk A, Davis S, De Moor B, Brazma A, Huber W
#' (2005). “BioMart and Bioconductor: a powerful link between biological
#' databases and microarray data analysis.” Bioinformatics, 21, 3439–3440.
#'
#' Wickham H, François R, Henry L, Müller K (2022). dplyr: A Grammar of Data
#' Manipulation.
#' \href{https://CRAN.R-project.org/package=dplyr}{Link}
#'
#' Müller K, Wickham H (2022). tibble: Simple Data Frames.
#' \href{https://CRAN.R-project.org/package=tibble}{Link}
#'
#' Lawrence M, Huber W, Pagès H, Aboyoun P, Carlson M, Gentleman R, Morgan M,
#' Carey V (2013). “Software for Computing and Annotating Genomic Ranges.” PLoS
#' Computational Biology, 9. doi: 10.1371/journal.pcbi.1003118,
#'
#' @export
#' @import utils
#' @importFrom stats na.omit
#' @importFrom tibble column_to_rownames
#' @importFrom pheatmap pheatmap
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

generateMatrix <- function(sfSeq = NA,
                           refTrp = NA,
                           sampleNames = NA,
                           type = "salmon",
                           keytype = "TXNAME",
                           element = "abundance",
                           species = NA,
                           release = NA,
                           outputCSV = FALSE,
                           nameCSV = "matrixCSV",
                           heatmap = T,
                           head = T) {

    if (typeof(sfSeq) != "character") {
        stop("Please provide valid .sf quantification files.")
    } else {
        ;
    }


    # Make a txdb object for the provided species.
    # If no reference transcriptome provided, obtain the annotation from Ensembl
    if (is.na(refTrp)) {

        # If no species name is provided, stop with an error message
        if (is.na(species)) {

            stop("No valid reference transcriptome file or species provided.")
        }

        if (is.na(release)) {
            refTrp <- obtainGTF(species = species, download = T)
        } else {
            refTrp <- obtainGTF(species = species, wantedVersion = release,
                                download = T)
        }

        # Unzip the gtf file
        system(paste("gunzip", refTrp))

        # Reformat the file name used for further process
        refTrp <- sub(".gz", "", refTrp)

    }

    txi <- txdbObj(sfSeq = sfSeq,
                   refTrp = refTrp,
                   species = species,
                   release = release,
                   key = keytype)

    # Create temporary file
    tmpFile <- "tmp.csv"

    if (element == "abundance") {
        write.csv(txi$abundance, file = tmpFile)
        }
    else if (element == "count") {
        write.csv(txi$count, file = tmpFile)
    }
    else if (element == "length") {
        write.csv(txi$length, file = tmpFile)
    } else {
        stop("Not an valid element. Please use abundance, count, or length.")
    }

    tmpMatrix <- read.csv(tmpFile)

    annot <- data.frame(tmpMatrix)

    colnames(annot) <- c("Gene ID", sampleNames)

    # A function used as a parameter of the "apply" function to filter out the
    # unexpressed genes.
    unEx <- function(x) {
        return(! all(x == 0))
    }

    # Get rid of unexpressed genes
    matrix <- annot[apply(annot, 1, unEx), ]

    if (outputCSV == TRUE) {
        write.csv(matrix, file = paste0(nameCSV, ".csv"))
    }

    # Delete the temporary file
    unlink("tmp.csv")
    unlink(refTrp)
    unlink(paste0(refTrp, ".gz"))

    # Format the matrix
    formatted <- tibble::column_to_rownames(matrix, "Gene ID")

    if (heatmap == T) {
        if (head == T) {

            (plot <- head(formatted))
            pheatmap::pheatmap(mat = plot,
                               display_numbers = T,
                               number_color = "black",
                               hclustfun = hclust)
        } else {
            plot <- formatted
        }

        pheatmap::pheatmap(mat = plot,
                           display_numbers = T,
                           number_color = "black",
                           hclustfun = hclust)
    }
    return(formatted)
}

# [END]
