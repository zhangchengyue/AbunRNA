library("GenomicFeatures")
library("tximport")

generateMatrix <- function (sfSeq, refSeq=NA, sampleNames, type="salmon", species=NA,
                            release=NA, abunCSV){

    # ... These arguments are the names of salmon output sf files.
    # sampleNames is charactor vector of sample names of the raw sequences.
    # type: Default tool for generating abundances is set to salmon.
    # species: Reference transcriptome of a species (the specific name)
    # release: The release of the reference. Default is the current release.
    # abunCSV: a character string for the name of the output abundance CSV

    # Make a txdb object for the provided species.

    if (is.na(refSeq)){
        if (is.na(species)){
            stop("No valid reference transcriptome file or species provided.")
        }
        txdb <- GenomicFeatures::makeTxDbFromEnsembl(species, release=release)
    }else{
        txdb <- GenomicFeatures::makeTxDbFromGFF(refSeq)
    }

    # Obtain txname keys
    k <- biomaRt::keys(txdb, keytypes = "TXNAME")
    txGene <- select(txdb, keys = k, columns = "TXNAME", keytype = "GENEID")
    txGene <- txGene[, c("TXNAME", "GENEID")]
    # return(txGene)

    txi <- tximport(sfSeq, type="salmon", tx2gene=txGene)
    abunFile <- paste(abunCSV, ".csv")
    write.csv(txi$abundance, file = abunFile)
    abundance <- read.csv(abunFile)
    abunAnnot <- data.frame(abundance)
    colnames(abunAnnot) <- c("Gene ID", sampleNames)

    # Get rid of unexpressed genes
    abunAnnot[apply(abunAnnot, 1, function(x){ !all(x==0)}),]
    return(abunAnnot)

}

# Try if this works.
# TODO: Make sure to modify quant.sf files to inst exdata subdirectory,
# and access them using system.file() in the test script.

sfSe <- c("cw1_quants/quant.sf", "cl1_quants/quant.sf", "cg1_quants/quant.sf")
names <- c("WT_WC_1", "lf_WC_1", "gf_WC_1")
output <- generateMatrix(sfSeq = sfSe,
                         refSeq="Caenorhabditis_elegans.WBcel235.107.gtf",
                         names, abunCSV = "abunOUT")
View(output)


