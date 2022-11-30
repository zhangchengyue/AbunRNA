#' Expression Abundance Matrix From C.elegans Quantification Result
#'
#' An RNAseq experiment conducted using C.elegans from 2021-2022 in Canada.
#'
#' @source Lunenfeld Tanenbaum Research Institute
#'
#' @format A data frame with 21390 rows and 3 variables:
#' \describe{
#'    \item{WT_WC_1}{ild type C.elegans sample sequence}
#'    \item{lf_WC_1}{loss of function C.elegans sample sequence}
#'    \item{gf_WC_1}{gain of function C.elegans sample sequence}
#' }
#'
#' @examples
#' \dontrun{
#' abundanceMatrix
#' }
"abundanceMatrix"

#' Genotypes for C.elegans samples
#'
#' An example file indicating the genotype conditions of C.elegans samples. Used
#' for displaying `plotPCA()` function utility.
#'
#' @source Lunenfeld Tanenbaum Research Institute
#'
#' @format A data frame with 3 rows and 1 variable:
#' \describe{
#'    \item{genotype}{indicating the genotype of the sample i.e. wild type,
#'    loss of function, gain of function}
#'    \item{samples}{sample names}
#'}
#'
#' @examples
#' \dontrun{
#' conditionsDF
#' }
"conditionsDF"

# [END]
