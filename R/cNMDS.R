#' Constrained Nonmetric Multidimensional Scaling
#'
#' @param formula Model formula where the left-hand side must be
#'     dissimilarities, and right-hand side the constraints.
#' @param data Data frame to find the terms on the right-hand side of
#'     the formula.
#' @param k Number of dimensions in NMDS.
#'
#' @importFrom stats cmdscale
#' @importFrom vegan monoMDS, MDSaddpoints
#'
#' @export
`cNMDS` <-
    function(formula, data, k = 2)
{
    ## step 1: constrained dissimilarities
    cdis <- distconstrain(formula, data)
    ## step 2: NMDS of constrained dissimilarities
    m0 <- cmdscale(cdis, k = k)
    cdis[] <- rank(cdis, ties = "min")
    sol <- monoMDS(cdis, m0, k = k)
    ## step 3: Constrained ommunity ordination
    dis <- eval(formula[[2]])
    dis[] <- rank(dis, ties = "min")
    m2 <- MDSaddpoints(sol, as.matrix(dis))
    ## construct result object
    sol$constraints <- sol$points
    sol$points <- m2$points
    sol$stress <- sol$stress + m2$deltastress
    class(sol) <- c("cNMDS", class(sol))
    sol
}
