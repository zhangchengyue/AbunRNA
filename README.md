
<!-- README.md is generated from README.Rmd. Please edit that file -->

# AbunRNA

An R package developed to focus on the commonly used transcript abuncance analysis pipelines including indexing, quantification, and to the downstream expression analysis.This R-based pipeline can access the available cDNA, DNA and transcriptome data from Ensembl.

<!-- badges: start -->
<!-- badges: end -->

## Installation

A paragraph that describes the purpose of your R package. Explain how
your package add to or improve a current work flow in bioinformatics or
computational biology (i.e., how is it unique?, what issue does it
address?). Finally, include the R version (not RStudio version) and
platform (Mac, Windows, Linux (Debian, Fedora/Redhat, Ubuntu)), used to
develop the package. There should be no `Shiny` implementation at this
point. You may obtain this information by running
`utils::sessionInfo()`. E.g., <br> <br> <br>

`AbunRNA` is an R package to demonstrate components of a simple R
package. This includes the main components: DESCRIPTION, NAMESPACE, man
subdirectory and R subdirectory. Additionally, LICENSE, README and
subdirectories vignettes, tests, data and inst are also explored. The
package is targeted for BCB410H (Applied Bioinformatics) students, who
are to define a useful tool for the analysis of biological data in the
format of a public R package housed on GitHub. The scope of the R
package is to add to or improve a current work flow in bioinformatics or
computational biology. The tool should contain functions to perform
analysis of biological data and to produce a compelling graphical
output, ideally to support for exploratory analysis. The
`TestingPackage` package was developed using
`R version 4.1.1 (2021-08-10)`,
`Platform: x86_64-apple-darwin17.0 (64-bit)` and
`Running under: macOS Big Sur 11.2`.

## Installation

To install the latest version of the package:

``` r
require("devtools")
install_github("zhangchengyue/AbunRNA", build_vignettes = TRUE)
library("AbunRNA")
```

To run the Shiny app:

``` r
runAbunRNA()
```

## Overview

Provide the following commands, customized to your R package. Then
provide an overview to briefly describe the main components of the
package. Include one image illustrating the overview of the package,
that shows the inputs and outputs. Ensure the image is deposited in the
correct location, as discussed in class. Point the user to vignettes for
a tutorial of your package. E.g., <br> <br> <br>

``` r
ls("package:TestingPackage")
data(package = "TestingPackage") 
browseVignettes("TestingPackage")
```

`TestingPackage` contains 4 functions to demonstrate components of a
simple R package. The *InfCriteriaCalculation* function calculates the
information criteria values. Specifically, Bayesian information
criterion (BIC), Akaike information criterion (AIC) and Integrated
Complete Likelihood (ICL) are calculated. The *InfCriteriaPlot*
generates a plot of information criteria values. *NormFactors* is a
function that calculates normalization factors via Trimmed Mean of
M-values (TMM). The *runTestingPackage* is the function that launches
the shiny app for this package. The package also contains two RNA
sequencing datasets, called GeneCounts and GeneCounts2. Refer to package
vignettes for more details. An overview of the package is illustrated
below.

![](./inst/figures/README-pressure-1.png)
