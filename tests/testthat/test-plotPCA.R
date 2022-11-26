library(AbunRNA)
library(testthat)
test_that("Check for input matrix", {

    expect_error(plotPCA(scaleIt = TRUE, conditions = conditionsDF,
                         col = "genotype"),
                 "Invalid input matrix.")
})

#> Test passed :)


testthat::test_that("Check for valid input arguments", {

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
    countMatrix <- generateMatrix(sfSeq = sfSe,
                                  species = "Caenorhabditis Elegans",
                                  sampleNames = samples,
                                  heatmap = FALSE,
                                  outputCSV = FALSE)


    expect_error(plotPCA(matrix = countMatrix, scaleIt = TRUE,
                         conditions = conditionsDF,
                         col = "genotype", x = 100, y = 100),
                 "Invalid x and y variables")

    expect_error(plotPCA(matrix = countMatrix, scaleIt = TRUE,
                         conditions = conditionsDF, x = 1, y = 2),
                 "col argument not provided")

    expect_error(plotPCA(matrix = countMatrix, scaleIt = TRUE,
                         conditions = conditionsDF,
                         col = "notexists", x = 1, y = 2),
                 "col argument not a column of conditions.")
})

#> Test passed :)

