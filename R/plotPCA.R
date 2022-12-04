#' Principle component analysis on input matrix
#'
#' A function that performs principle component analysis on input matrix, and
#' plot the graph indicaiting the results, grouped in category given by
#' the user.
#'
#' @param matrix A count matrix. Default is NULL.
#' @param scaleIt A boolean indicating whether to scale the variables (divide by
#'     standard deviation). Default value is set to TRUE.
#' @param conditions A data frame indicating the conditions of the sample.
#'     Default is NULL.
#' @param col A string indicating the category where observations should be
#'     grouped into. This should be a column name in "conditions" data frame.
#'     Default is NA.
#' @param x An integer indicating the x variable (e.g. x = 1 would extrace PC1
#'     from principle component table). Default would be PC1.
#' @param y An integer indicating the y variable. Default would be PC2.
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
#' library("ggfortify") # Load ggfortify for plotting PCA
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
#' @importFrom stats prcomp
#' @importFrom stats hclust
#' @importFrom ggplot2 autoplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @import ggfortify

plotPCA <- function(matrix = NULL, scaleIt = TRUE,
                    conditions = NULL, col = NA, x = 1, y = 2) {

    if (is.null(matrix)) {
        stop("Invalid input matrix.")
    } else {
        ;
    }

    if (x == y){
        stop("x and y variables cannot be identical.")
    } else {
        ;
    }

    if (is.na(col) && (!is.null(conditions))) {
        stop("col argument not provided")
    } else {
        ;
    }

    if (! (col %in% colnames(conditions))) {
        stop("col argument not a column of conditions.")
    } else {
        ;
    }

    transpose <- as.data.frame(t(matrix))
    rownames(transpose) <- NULL

    # Filter
    transpose <- transpose[, colSums(transpose) > 0]

    pca <- stats::prcomp(x = transpose, scale = scaleIt)

    pcaR <- data.frame(pca$rotation)
    if (! x %in% seq(along = pcaR) || (! y %in% seq(along = pcaR))) {
        stop("Invalid x and y variables")
    } else {
        ;
    }


    if (is.null(conditions)) {
        cat("Conditions not provided. Data would not be categorized.")
        plot <- ggplot2::autoplot(pca, x = x, y = y)

    } else {

        groups <- conditions[ , colnames(conditions) == col]
        transpose$Conditions <- as.factor(groups)

        plot <- ggplot2::autoplot(pca, data = transpose,x = x, y = y) +
            ggplot2::geom_point(ggplot2::aes(color = groups))
    }

    return(list("PCA" = pcaR, "Plot" = plot))
}

# [END]
