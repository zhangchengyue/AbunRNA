library(AbunRNA)

test_that("Quantification files not valid", {
    expect_error(txdbObj(sfSeq = 123,
                                refTrp = NA,
                                key = "TXNAME",
                                species = NA,
                                release = NA),
                 "Please provide valid .sf quantification files.")

    expect_error(txdbObj(sfSeq = TRUE,
                         refTrp = NA,
                         key = "TXNAME",
                         species = NA,
                         release = NA),
                 "Please provide valid .sf quantification files.")

    expect_error(txdbObj(sfSeq = list(1, 2, 3),
                         refTrp = NA,
                         key = "TXNAME",
                         species = NA,
                         release = NA),
                 "Please provide valid .sf quantification files.")

    expect_error(txdbObj(sfSeq = NA,
                         refTrp = NA,
                         key = "TXNAME",
                         species = NA,
                         release = NA),
                 "Please provide valid .sf quantification files.")
})

#> Test pass :)

test_that("Make sure key is in keytypes", {

    annot <- obtainGTF(species = "Caenorhabditis elegans",
                        download = T)

    system(paste("gunzip", annot))

    # Reformat the file name used for further process
    annot <- sub(".gz", "", annot)

    txdb <- GenomicFeatures::makeTxDbFromGFF(annot)

    # Obtain keytypes
    allKeys <- biomaRt::keytypes(txdb)
    err <- paste("Key doesn't exist. Available keytypes: ",
                 allKeys, sep = ",", collapse = ", ")


    cw1Quants <- system.file("extdata",
                             "cw1_quants",
                             "quant.sf",
                             package = "AbunRNA")

    cl1Quants <- system.file("extdata",
                             "cl1_quants",
                             "quant.sf",
                             package = "AbunRNA")

    cg1Quants <- system.file("extdata",
                             "cg1_quants",
                             "quant.sf",
                             package = "AbunRNA")

    sfSe <- c(cw1Quants, cl1Quants, cg1Quants)

    expect_error(txdbObj(sfSeq = sfSe,
                         refTrp = annot,
                         key = "NOTEXIST",
                         species = NA,
                         release = NA), err)
    unlink(annot)
    unlink(paste0(annot, ".gz"))

})

#> Test pass :)
