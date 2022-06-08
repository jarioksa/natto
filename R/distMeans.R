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
#' @param d Distances as a \code{dist} object
#' @return Distances to all other points from a point that is in the
#'     centre of the coordinates generating the distances.
#' @export
`distMeans` <-
    function(d)
{
    x <- as.matrix(d^2/2)
    ## Gower double centring
    x <- sweep(x, 2L, colMeans(x), check.margin = FALSE)
    x <- sweep(x, 1L, rowMeans(x), check.margin = FALSE)
    ## d is now Gower double-centred, and we need backtransform points
    ## at zero back to distances to all other points. For full matrix
    ## this would be sqrt(2*d - outer(diag(d), diag(d), "+")), but we
    ## only need the centroid for a zero-row (d == 0).
    sqrt(diag(-x))
}

