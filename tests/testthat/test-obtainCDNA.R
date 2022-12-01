library(AbunRNA)
library(testthat)
test_that("Extracts the correct cDNA gz file for C.elegans, version 107", {
    file <- obtainCDNA(species = "Caenorhabditis Elegans",
                       wantedVersion=107,
                       download = F)
    expect_equal(file, "Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz")
})

#> Test passed :)

test_that("Raise error if species name is not provided", {
    expect_error(obtainCDNA(wantedVersion=107, download = F),
                 "Species name not provided.")
})

#> Test passed :)

test_that("Raise error if species name is incorrect or doesn't exist", {
    expect_error(obtainCDNA(species = "Nonexist Species", wantedVersion=107),
                 "Species name cannot be recognized. Please try again.")
})

#> Test passed :)



