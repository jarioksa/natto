### Centroid of distances cannot be found directly as the mean of
### distances, but we need to use Gower double centralized data where
### the real mean of points is 0, and then backtransform to
### distances. This is based on Euclidean geometry -- but so is all
### Hmsc.
#' Mean of Distances
#'
#' Mean of distances is defined as the distance of each point to the
#' mean of coordinates generating the distance.
#'
#' @details Function is analagous to \code{\link{colMeans}} or
#'     \code{\link{rowMeans}} and returns values that are at the mean
#'     of distances of each row or column of a symmetric distance
#'     matrix. Such means cannot be directly found as marginal means
#'     of distance matrix, but they must be found after Gower double
#'     centring (Gower 1966). After double centring, the means are
#'     zero, and when backtransformed to original distances, these
#'     give the mean distances. When added to the original distances,
#'     the metric properties are preserved. For instance, adding
#'     centres to distances will not influence results of metric
#'     scaling, or the rank of spatial Euclidean distances. The method
#'     is based on Euclidean geometry, but also works for
#'     non-Euclidean dissimilarities. However, the means of very
#'     strongly non-Euclidean indices may be imaginary, and given as
#'     \code{NaN}.
#'
#' @references Gower, J.C. (1966) Some distance properties of latent
#'     root and vector methods used in multivariate analysis.
#'     \emph{Biometrika} \bold{53}, 325-328.
#'
#' @author Jari Oksanen.
#'
#' @param d Distances as a \code{dist} object
#' @param addcentre Add distances to the centroid as the first item in
#'     the distance matrix. If \code{FALSE} only return mean
#'     distances.
#' @param label Label for the centroid when \code{addcentre = TRUE}.
#' @param formula,data Formula where the left-hand-side is the
#'     dissimilarity structure, and right-hand-side defines the mean
#'     from which the dissimilarities are calculated. The terms in the
#'     right-hand-side can be given in \code{data}.
#' @param \dots Other parameters (ignored).
#'
#' @return Either distances to all other points from a point that is
#'     in the centroid of the coordinates generating the distances, or
#'     the input dissimilarity matrix where the mean distances are
#'     added as the first observation.
#' @examples
#' ## Euclidean distances to the mean of coordinates ...
#' xy <- matrix(runif(5*2), 5, 2)
#' dist(rbind(xy, "mean" = colMeans(xy)))
#' ## ... are equal to distMeans ...
#' distMeans(dist(xy))
#' ## ... but different from mean of distances
#' colMeans(as.matrix(dist(xy)))
#' ## adding mean distance does not influence PCoA of non-Euclidean
#' ## distances (or other metric properties)
#' data(spurn)
#' d <- canneddist(spurn, "bray")
#' m0 <- cmdscale(d, eig = TRUE)
#' mcent <- cmdscale(distMeans(d, addcentre=TRUE), eig = TRUE)
#' ## same non-zero eigenvalues
#' zapsmall(m0$eig)
#' zapsmall(mcent$eig)
#' ## distMeans are at the origin of ordination
#' head(mcent$points)
#' @export
`distMeans` <-
    function(...)
        UseMethod("distMeans")

#' @rdname distMeans
#' @export
`distMeans.default` <-
    function(d, addcentre = FALSE, label = "centroid", ...)
{
    x <- as.matrix(d^2/2)
    ## Gower double centring
    x <- sweep(x, 2L, colMeans(x), check.margin = FALSE)
    x <- sweep(x, 1L, rowMeans(x), check.margin = FALSE)
    ## x is now Gower double-centred, and we need backtransform points
    ## at zero back to distances to all other points. For full matrix
    ## this would be sqrt(2*x - outer(diag(x), diag(x), "+")), but
    ## after double-centring the centre (or the mean) is zero.
    cnt <- sqrt(diag(-x))
    if (addcentre) {
        att <- attributes(d)
        cnt <- c(cnt, d)
        att$Size <- att$Size + 1L
        if (!is.null(att$Labels))
            att$Labels <- c(label, att$Labels)
        attributes(cnt) <- att
    }
    cnt
}

#' @rdname distMeans
#' @export
`distMeans.formula` <-
    function(formula, data, ...)
{
    x <- vegan:::initDBRDA(eval(formula[[2]]))
    mm <- model.matrix(delete.response(terms(formula)), data)
    sqrt(diag(vegan:::ordPartial(x, mm)$Y))
}
