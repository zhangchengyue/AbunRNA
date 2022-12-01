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
#'    \item{genotype}{Genotypes of each sample sequence. "wt" is wild type,
#'    "lf" is loss of function, "gf" is gain of function}
#'    \item{samples}{Sample sequences}
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
#'    \item{cg1}{gain of function sequence sample}
#'    \item{cg2}{gain of function sequence sample}
#'    \item{cg3}{gain of function sequence sample}
#'    \item{sg1}{gain of function sequence sample}
#'    \item{sg2}{gain of function sequence sample}
#'    \item{sg3}{gain of function sequence sample}
#'    \item{cl1}{loss of function sequence sample}
#'    \item{cl2}{loss of function sequence sample}
#'    \item{cl3}{loss of function sequence sample}
#'    \item{sl1}{loss of function sequence sample}
#'    \item{sl2}{loss of function sequence sample}
#'    \item{sl3}{loss of function sequence sample}
#'    \item{cw1}{wild type sequence sample}
#'    \item{cw2}{wild type sequence sample}
#'    \item{cw3}{wild type sequence sample}
#'    \item{sw1}{wild type sequence sample}
#'    \item{sw2}{wild type sequence sample}
#'    \item{sw3}{wild type sequence sample}
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
