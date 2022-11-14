#' Principle component analysis on input matrix
#'
#'
#' A function that performs principle component analysis on input matrix, and
#' plot the graph indicaiting the results, grouped in category given by
#' the user.
#'
#' @param mat A count matrix
#' @param scale A boolean indicating whether to scale the variables (divide by
#'     standard deviation). Default value is set to TRUE
#' @param conditions A data frame indicating the conditions of the sample
#' @param col A string indicating the category where observations should be
#'     grouped into. This should be a column name in "conditions" data frame.
#' @param x An integer indicating the x variable (e.g. x = 1 would extrace PC1
#'     from principle component table)
#' @param y An integer indicating the y variable
#'
#'
#' @return Returns a principle component analysis table, and a PCA plot for
#'      x and y variables
#'
#' @examples
#' cw1_quants <- system.file("extdata",
#'                           "cw1_quants",
#'                           "quant.sf",
#'                           package = "AbunRNA")
#' cl1_quants <- system.file("extdata",
#'                           "cl1_quants",
#'                           "quant.sf",
#'                           package = "AbunRNA")
#' cg1_quants <- system.file("extdata",
#'                           "cg1_quants",
#'                           "quant.sf",
#'                           package = "AbunRNA")
#' sfSe <- c(cw1_quants, cl1_quants, cg1_quants)
#'
#' # Name the samples correspondingly
#' samples <- c("WT_WC_1", "lf_WC_1", "gf_WC_1")
#' matrix <- generateMatrix(sfSeq = sfSe,
#'                       species = "Caenorhabditis elegans",
#'                       release = 107,
#'                       sampleNames = samples,
#'                       outputCSV = TRUE,
#'                       abunCSV = "abunOUT")
#' conditionsDF <- read.csv(system.file("extdata", "conditions.csv")
#'
#' # Note that the col argument for plotPCA should be a column name
#' # in conditionsDF
#' View(conditionsDF)
#'
#' pca <- plotPCA(matrix, scale = TRUE,
#'                        conditions = conditionsDF,
#'                        col = "genoptype")
#'
#' @references
#' Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis.
#' Springer-Verlag New York. ISBN 978-3-319-24277-4
#' \href{https://ggplot2.tidyverse.org}
#'
#' R Core Team (2013). R: A language and environment for statistical computing.
#' R Foundation for Statistical Computing, Vienna, Austria. ISBN 3-900051-07-0
#' \href{http://www.R-project.org/}
#'
#' @export
#' @import ggplot2
#' @import stats

plotPCA <- function(mat, scale = TRUE, conditions, col, x = PC1, y = PC2){
    # Principle component of abuncance data
    pca <- stats::prcomp(mat, scale = scale)

    pcaR<-data.frame(pca$rotation)

    for (i in seq(length(colnames(pcaR)))){
        assign(paste0("PC", i), pcaR[, i])
    }

    Groups <- conditions[, colnames(conditions) == col]
    PCx <- x
    PCy <- y

    ggplot2::ggplot(pcaR,
                    ggplot2::aes(x=PCx, y=PCy, fill = Groups)) +
        ggplot2::geom_point(shape = 21, col = "black")
    return(pcaR)
}

# [END]
