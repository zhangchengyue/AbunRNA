#' Genotypes for C.elegans larger samples
#'
#' A condition matrix for samples in abunMatrix, describing the genotype
#' of each sample.
#' 18 samples are included.
#'
#' @source Lunenfeld Tanenbaum Research Institute
#'
#' @format A data frame with 18 rows and 2 variables.
#' \describe{
#'    \item{lf}{gain of function genotype}
#'    \item{cl_x/sl_x}{loss of function genotype}
#'    \item{cw_x/sw_x}{wild type genotype}
#' }
#'
#' @examples
#' \dontrun{
#' bigCond
#' }
"bigCond"

#' Larger Sample Size Expression Abundance Matrix From C.elegans Quantification Result
#'
#' An RNAseq experiment conducted using C.elegans from 2021-2022 in Canada.
#' 18 samples are included.
#'
#' @source Lunenfeld Tanenbaum Research Institute
#'
#' @format A data frame with 21390 rows and 18 variables.
#' \describe{
#'    \item{cg_x/sg_x}{gain of function sequences}
#'    \item{cl_x/sl_x}{loss of function sequences}
#'    \item{cw_x/sw_x}{wild type sequences}
#' }
#'
#' @examples
#' \dontrun{
#' abunMatrix
#' }
"abunMatrix"

#' Small Sample Size Expression Abundance Matrix From C.elegans Quantification Result
#'
#' An RNAseq experiment conducted using C.elegans from 2021-2022 in Canada.
#' Only 3 samples are included.
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
#' countMatrix
#' }
"countMatrix"

#' Genotypes for C.elegans small samples
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
