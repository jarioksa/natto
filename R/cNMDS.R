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
#'     \code{\link[vegan]{metaMDS}}.
#'   \item Add points to constrained ordination using unconstrained
#'     community dissimilarities with \code{\link[vegan]{MDSaddpoints}}.
#' } % end enumerate
#'
#' The algorithms creates analogous result to constrained eigenvector
#' ordination such as distance-based RDA (\code{\link[vegan]{dbrda}}).
#' The constrained NMDS is similar to linear combination (LC) scores
#' and the site NMDS analogous to to weighted averages (WA) scores of
#' db-RDA. The \code{cNMDS} site scores are found by adding the sites
#' with their community dissimilarities as new points to constrained
#' ordination. Instead of raw dissimilarities, both the constrained
#' ordination and the site ordination ("new points") are based on
#' ranks of dissimilarities (with ties). The ranks give the same
#' result as raw dissimilarities within one set of dissimilarities
#' (constrained, community), but with them these constrained and
#' community ordination are rank-order similar.
#'
#' Although the idea of the algorithm is simple and obvious, the
#' results are often far from satisfactory. Quite often the community
#' dissimilarities and constrained dissimilarities differ from each
#' other so strongly that the new added points (community) are
#' expelled from the constrained configuration and are located at the
#' outskirts of the ordination instead of being mixed with the
#' constrained points.
#'
#' @return Function returns an object of class \code{"cNMDS"} that
#'     inherits from \pkg{vegan} functions
#'     \code{\link[vegan]{metaMDS}} and
#'     \code{\link[vegan]{monoMDS}}. but it has rudimentary
#'     \code{scores} and \code{plot} methods. The \code{plot} only
#'     shows one set of scores, and it is best to build graphics with
#'     pipes (see Examples). The result object also contains biplot
#'     arrows for econtinuous constraints and centroids for factor
#'     constraints. These are based on community ordination (sites)
#'     instead of constraints, and they are added with
#'     \code{\link[vegan]{envfit}}.
#'
#' @seealso \code{\link{distconstrain}}, \code{\link[vegan]{dbrda}},
#'     \code{\link[vegan]{monoMDS}}, \code{\link[vegan]{metaMDS}}.
#'
#' @section Warning:
#'
#' The function is provided as an example of the
#' algorithm. The results are often poor, and the function should not
#' be used for any other purposes than inspecting the method.
#'
#' @examples
#' if (require(vegan)) {
#' data(mite, mite.env)
#' dis <- canneddist(mite, "chord")
#' mod <- cNMDS(dis ~ WatrCont + SubsDens + Topo + Shrub, mite.env)
#' mod
#' plot(mod, type = "p") |>
#'    points("constraints", pch = 16, col = 2) |>
#'    text("biplot", col = 4) |>
#'    text("centroids", col = 4)
#' }
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
#' @importFrom vegan metaMDS MDSaddpoints envfit scores
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
    cdis[] <- rank(round(cdis, 12), ties.method = "min")
    sol <- metaMDS(cdis, m0, k = k, trace = FALSE)
    ## step 3: Constrained community ordination
    dis <- eval(formula[[2]])
    dis[] <- rank(round(dis, 12), ties.method = "min")
    m2 <- MDSaddpoints(sol, as.matrix(dis))
    ## This was the last step: the rest is janitorial and adding candies
    terms <- delete.response(terms(formula, data = data))
    mf <- model.frame(terms, data = data)
    ef <- envfit(m2, mf)
    if (!is.null(ef$vectors)) {
        sol$biplot <- scores(ef, "vectors")
        attr(sol$biplot, "score") <- "biplot"
    }
    if (!is.null(ef$factors))
        sol$centroids <- scores(ef, "factors")
    ## construct result object: inherits from vegan::monoMDS
    sol$constraints <- sol$points
    sol$sites <- m2$points
    sol$stress <- sol$stress + m2$deltastress
    attr(sol$points, "pc") <- FALSE
    sol$call <- match.call()
    sol$model <- "Constrained non-metric"
    sol$distmethod <- attr(dis, "method")
    sol$distcall <- NULL
    sol$iters <- m2$iters
    sol$icause <- m2$cause
    sol$iscal <- FALSE
    class(sol) <- c("cNMDS", "monoMDS")
    sol
}

#' @rdname cNMDS
#' @param x \code{cNMDS} result object.
#' @param display Kind of scores to display. Can be one or several of
#'     \code{"sites"}, \code{"constraints"}, \code{"biplot"},
#'     \code{"centroids"}, or alternative \code{"all"} for all these.
#' @export
`scores.cNMDS` <-
    function(x, display = "sites", ...)
{
    scores <- c("sites", "constraints", "biplot", "centroids")
    display <- match.arg(display, c("all", scores), several.ok = TRUE)
    if ("all" %in% display)
        display <- scores
    keep <- display %in% names(x)
    display <- display[keep]
    out <- x[display]
    if (is.list(out) && length(out) == 1)
        out <- out[[1]]
    out
}

#' @rdname cNMDS
#' @param type Either \code{"t"}ext, \code{"p"}oints or \code{"n"}one.
#' @param \dots Other arguments passed to graphics functions.
#'
#' @export
`plot.cNMDS` <-
    function(x, display = "sites", type = "p", ...)
{
    if (length(display) > 1)
        stop("only one item can be used in plot: use text(), points() in pipe for others")
    out <- scores(x, display = "all")
    xlim <- range(sapply(out, function(z) z[,1]))
    ylim <- range(sapply(out, function(z) z[,2]))
    plt <- scores(x, display = display)
    suppressMessages(ordiplot(plt, type = type, xlim = xlim, ylim = ylim, ...))
    class(out) <- c("cNMDS", "ordiplot")
    invisible(out)
}
