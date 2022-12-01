library(AbunRNA)
library(testthat)
test_that("Quantification files not valid", {
    expect_error(generateMatrix(sfSeq = 123,
                                refTrp = NA,
                                sampleNames = "Caenorhabditis Elegans",
                                type = "salmon",
                                keytype = "TXNAME",
                                species = NA,
                                release = NA,
                                outputCSV = FALSE,
                                abunCSV = "abunCSV",
                                heatmap = T,
                                head = T),
                 "Please provide valid .sf quantification files.")

    expect_error(generateMatrix(sfSeq = TRUE,
                                refTrp = NA,
                                sampleNames = "Caenorhabditis Elegans",
                                type = "salmon",
                                keytype = "TXNAME",
                                species = NA,
                                release = NA,
                                outputCSV = FALSE,
                                abunCSV = "abunCSV",
                                heatmap = T,
                                head = T),
                 "Please provide valid .sf quantification files.")

    expect_error(generateMatrix(sfSeq = list(1, 2, 3),
                                refTrp = NA,
                                sampleNames = "Caenorhabditis Elegans",
                                type = "salmon",
                                keytype = "TXNAME",
                                species = NA,
                                release = NA,
                                outputCSV = FALSE,
                                abunCSV = "abunCSV",
                                heatmap = T,
                                head = T),
                 "Please provide valid .sf quantification files.")

    expect_error(generateMatrix(refTrp = NA,
                                sampleNames = "Caenorhabditis Elegans",
                                type = "salmon",
                                keytype = "TXNAME",
                                species = NA,
                                release = NA,
                                outputCSV = FALSE,
                                abunCSV = "abunCSV",
                                heatmap = T,
                                head = T),
                 "Please provide valid .sf quantification files.")
})

#> Test pass :)

test_that("Generates correct count matrix", {
    expectedMatrix <- list(WT_WC_1 = c(7.808453, 5.207178, 5.541982,
                                       3.897999, 4.190805, 3.571827),
                           lf_WC_1 = c(10.264917, 8.038582, 13.091974,
                                       6.805610, 5.422508, 6.277508),
                           gf_WC_1 = c(8.107014, 6.298857, 7.474552, 5.752922,
                                       3.959048, 6.189619))

    expectedMatrix <- data.frame(expectedMatrix)

    rownames(expectedMatrix) <- c("WBGene00000001", "WBGene00000002",
                                  "WBGene00000003", "WBGene00000004",
                                  "WBGene00000005", "WBGene00000006")

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
    samples <- c("WT_WC_1", "lf_WC_1", "gf_WC_1")
    actualMatrix <- generateMatrix(sfSeq = sfSe,
                                   species = "Caenorhabditis Elegans",
                                   sampleNames = samples,
                                   outputCSV = TRUE)
    actualMatrix <- head(actualMatrix)
    expect_equal(expectedMatrix, actualMatrix)
})

#> Test pass :)

