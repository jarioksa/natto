#' Constrained Nonmetric Multidimensional Scaling
#'
#' Adds points based on community dissimilarities on the NMDS
#' based on constrained dissimilarities.
#'
#' The steps of algorithm are:
#' \enumerate{
#'   \item Find constrained dissimilarities for given constraints
#'     with \code{\link{distconstrain}}.
#'   \item Find NMDS of constrained dissimilarities with
#'     \code{\link[vegan]{monoMDS}} using metric scaling as initial
#'     configuration.
#'   \item Add points to constrained ordination using unconstrained
#'     community dissimilarities with \code{\link[vegan]{MDSaddpoints}}.
#' } % end enumerate
#'
#' @param formula Model formula where the left-hand side must be
#'     dissimilarities, and right-hand side the constraints.
#' @param data Data frame to find the terms on the right-hand side of
#'     the formula.
#' @param k Number of dimensions in NMDS.
#' @param add Type of additive constant used to avoid negative squared
#'     distances in \code{\link{distconstrain}}. Either
#'     \code{"lingoes"} or \code{"cailliez"} or \code{FALSE}.
#'
#' @importFrom stats cmdscale delete.response terms model.frame
#' @importFrom vegan monoMDS MDSaddpoints envfit scores
#'
#' @export
`cNMDS` <-
    function(formula, data, k = 2, add = FALSE)
{
    ## step 1: constrained dissimilarities
    cdis <- distconstrain(formula, data, add = add, squared = TRUE)
    if (any(cdis < 0)) {
        message("negative squared dissimilarities changed to 0")
        cdis[cdis < 0] <- 0
    }
    ## step 2: NMDS of constrained dissimilarities
    m0 <- cmdscale(sqrt(cdis), k = k)
    cdis[] <- rank(cdis, ties.method = "min")
    sol <- monoMDS(cdis, m0, k = k)
    ## step 3: Constrained community ordination
    dis <- eval(formula[[2]])
    dis[] <- rank(dis, ties.method = "min")
    m2 <- MDSaddpoints(sol, as.matrix(dis))
    ## This was the last step: the rest is janitorial and adding candies
    terms <- delete.response(terms(formula, data = data))
    mf <- model.frame(terms, data = data)
    ef <- envfit(m2, mf)
    if (!is.null(ef$vectors))
        sol$biplot <- scores(ef, "vectors")
    if (!is.null(ef$factors))
        sol$centroids <- scores(ef, "factors")
    ## construct result object: inherits from vegan::monoMDS
    sol$constraints <- sol$points
    sol$points <- m2$points
    sol$stress <- sol$stress + m2$deltastress
    attr(sol$points, "pc") <- FALSE
    sol$call <- match.call()
    sol$model <- "Constrained non-metric"
    sol$distmethod <- attr(dis, "method")
    sol$distcall <- NULL
    sol$iters <- m2$iters
    sol$icause <- m2$cause
    sol$iscal <- FALSE
    class(sol) <- c("cNMDS", class(sol))
    sol
}
