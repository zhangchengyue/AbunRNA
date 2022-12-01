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
#'    \item{cg_1}{gain of function sequence sample}
#'    \item{cg_2}{gain of function sequence sample}
#'    \item{cg_3}{gain of function sequence sample}
#'    \item{sg_1}{gain of function sequence sample}
#'    \item{sg_2}{gain of function sequence sample}
#'    \item{sg_3}{gain of function sequence sample}
#'    \item{cl_1}{loss of function sequence sample}
#'    \item{cl_2}{loss of function sequence sample}
#'    \item{cl_3}{loss of function sequence sample}
#'    \item{sl_1}{loss of function sequence sample}
#'    \item{sl_2}{loss of function sequence sample}
#'    \item{sl_3}{loss of function sequence sample}
#'    \item{cw_1}{wild type sequence sample}
#'    \item{cw_2}{wild type sequence sample}
#'    \item{cw_3}{wild type sequence sample}
#'    \item{sw_1}{wild type sequence sample}
#'    \item{sw_2}{wild type sequence sample}
#'    \item{sw_3}{wild type sequence sample}
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
