#' Principle component analysis on input matrix
#'
#' A function that performs principle component analysis on input matrix, and
#' plot the graph indicaiting the results, grouped in category given by
#' the user.
#'
#' @param matrix A count matrix
#' @param scalIt A boolean indicating whether to scale the variables (divide by
#'     standard deviation). Default value is set to TRUE
#' @param conditions A data frame indicating the conditions of the sample
#' @param col A string indicating the category where observations should be
#'     grouped into. This should be a column name in "conditions" data frame.
#' @param x An integer indicating the x variable (e.g. x = 1 would extrace PC1
#'     from principle component table)
#' @param y An integer indicating the y variable
#'
#' @return Returns a list containing the principle component analysis results,
#'     and a PCA plot for x and y variables.
#'
#' @examples
#' cw1Quants <- system.file("extdata",
#'                           "cw1_quants",
#'                           "quant.sf",
#'                           package = "AbunRNA")
#' cl1Quants <- system.file("extdata",
#'                           "cl1_quants",
#'                           "quant.sf",
#'                           package = "AbunRNA")
#' cg1Quants <- system.file("extdata",
#'                           "cg1_quants",
#'                           "quant.sf",
#'                           package = "AbunRNA")
#' sfSe <- c(cw1Quants, cl1Quants, cg1Quants)
#'
#' # Name the samples correspondingly
#' samples <- c("WT_WC_1", "lf_WC_1", "gf_WC_1")
#' countMatrix <- generateMatrix(sfSeq = sfSe,
#'                       species = "Caenorhabditis elegans",
#'                       release = 107,
#'                       sampleNames = samples,
#'                       outputCSV = FALSE)
#' graphPlot <- plotPCA(matrix = countMatrix, scaleIt = TRUE,
#'                      conditions = conditionsDF,
#'                      col = "genotype")
#' graphPlot$PCA
#' graphPlot$Plot
#'
#' @references
#' Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis.
#' Springer-Verlag New York. ISBN 978-3-319-24277-4
#' \href{ https://CRAN.R-project.org/package=ggplot2}{Link}
#'
#' R Core Team (2013). R: A language and environment for statistical computing.
#' R Foundation for Statistical Computing, Vienna, Austria. ISBN 3-900051-07-0
#' \href{http://www.R-project.org/}{Link}
#'
#' @export
#' @import ggplot2
#' @import stats

plotPCA <- function(matrix = NULL, scaleIt = TRUE,
                    conditions = NULL, col = NA, x = 1, y = 2) {

    if (is.null(matrix)) {
        stop("Invalid input matrix.")
    } else {
        ;
    }

    pca <- stats::prcomp(x = matrix, scale = scaleIt)

    pcaR <- data.frame(pca$rotation)
    if (! x %in% seq(along = pcaR) || (! y %in% seq(along = pcaR))) {
        stop("Invalid x and y variables")
    } else {
        ;
    }

    PCx <- pcaR[ , x]
    PCy <- pcaR[ , y]

    if (is.null(conditions)) {
        cat("Conditions not provided. Data would not be categorized.")
        plot <- ggplot2::ggplot(pcaR,
                                ggplot2::aes(x = PCx, y = PCy)) +
            ggplot2::geom_point(size = 5, shape = 21, col = "black")
    } else {
        if (is.na(col)) {
            stop("col argument not provided")
        } else {
            ;
        }

        if (! (col %in% colnames(conditions))) {
            stop("col argument not a column of conditions.")
        } else {
            ;
        }

        Groups <- conditions[ , colnames(conditions) == col]

        plot <- ggplot2::ggplot(pcaR,
                                ggplot2::aes(x = PCx, y = PCy, fill = Groups)) +
            ggplot2::geom_point(size = 5, shape = 21, col = "black")
    }

    return(list("PCA" = pcaR, "Plot" = plot))
}

# [END]
