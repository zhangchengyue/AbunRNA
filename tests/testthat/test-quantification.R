library(AbunRNA)
library(testthat)
test_that("Parse error message if fastq parameter is NA", {
    expect_error(quantification(species = "Caenorhabditis Elegans",
                  release = NA,
                  indexName = "celegansINDEX",
                  fastq = NA,
                  quantOut = "celegansQUANT"), "Invalid FASTQ files.")
})

#> Test pass :)
