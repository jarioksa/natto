#'
#' Convex Hull Enclosing a Given Proportion of Points
#'
#' Functions find a small convex hull or ellipse enclosing a given
#' proportion of points. The functions work by removing points from
#' the hull or ellipse one by one. The points are selected either so
#' that the area of the remaining hull is reduced as much as possible
#' at every step (de Smith, Goodrich & Longley 2007), or removing the
#' point that has the maximal total distance to all remaining
#' points or longest Mahalanobis distance to the group.

#' @encoding UTF-8
#'
#' @param pts Coordinates of points, a two-column matrix
#' @param keep Proportion of points kept
#' @param criterion Criterion to remove a point on the
#'     hull. \code{"area"} removes the point that reduces the area
#'     most, \code{"distance"} remove the point at the hull that has
#'     the largest total distance to all points on and within the
#'     hull, and \code{"mahalanobis"} remove the point with longest
#'     Mahalanobis distance to the centroid.
#'
#' @details
#'
#' Reduction of area is the only criterion that really is based the
#' area. The other methods are based on the distances to all points
#' within and on the hull or on the Mahalanobis distance. The
#' functions only work in 2D.
#'
#' The algorithms are naïve, and stepwise removal of single points
#' does not guarantee smallest possible final hull. Two outlier points
#' close to each other can protect each other from removal, but
#' removing both together would give a large reduction of the
#' hull. The \code{"distance"} criterion produces circular hulls, but
#' \code{"mahalanobis"} will better preserve the original elongation
#' of the configuration, although it rarely gives smaller areas than
#' \code{"area"} criterion.

#' @author Jari Oksanen
#'
#' @references de Smith, M.J., Goodchild, M.F. & Longley,
#'     P.A. (2007). \emph{Geospatial analysis: A comprehensive guide
#'     to principles, techniques and software tools}. Matador.
#'
#' @return A two-column matrix of coordinates of peeled hull or
#'     ellipse defining that can be plotted with
#'     \code{\link[graphics]{polygon}}, with attributes
#'     \code{"centre"} and \code{"area"}.
#'
#' @examples
#'
#' x <- matrix(c(rnorm(100), rnorm(100, sd=2)), nrow=100, ncol=2)
#' op <- par(mar=c(0,0,1,0), mfrow=c(2,2))
#' ## show first ten dropped points and hull keeping 50% of points
#' for (crit in c("area", "distance", "mahalanobis")) {
#'    plot(x, asp = 1, main = crit, xlab="", ylab="", axes=FALSE)
#'    for (p in c(100:90, 50)/100)
#'        polygon(peelhull(x, keep=p, crit=crit),
#'                col = adjustcolor("blue", alpha.f=0.04))
#' }
#' plot(x, asp = 1, main = "ellipse", xlab="", ylab="", axes=FALSE)
#' for (p in c(100:90, 50)/100)
#'    polygon(peelellipse(x, keep=p), col=adjustcolor("blue", alpha.f=0.04))
#'
#' par(op) # original graphical parameters
#'
#' @importFrom grDevices chull
#' @importFrom stats mahalanobis
#' @export
`peelhull` <-
    function(pts, keep = 0.9, criterion = c("area", "distance", "mahalanobis"))
{
    criterion <- match.arg(criterion)
    stopifnot(ncol(pts) == 2, keep <= 1, keep > 0)
    ndrop <- as.integer(nrow(pts) * (1-keep))
    ## remove one point on the hull, either to maximally reduce the
    ## area of the new hull, or a point that is most distant to all
    ## points on and within the hull.
    for (k in seq_len(ndrop)) {
        hull <- chull(pts)
        crit <- numeric(length(hull))
        for (i in seq_along(hull)) {
            crit[i] <- switch(criterion,
               "area" = polyarea(pts[-hull[i],][chull(pts[-hull[i],]),]),
               "distance" = -sum(sweep(pts, 2, pts[hull[i],])^2),
               "mahalanobis" = -mahalanobis(pts[hull[i],], colMeans(pts),
                                cov(pts))
               )
        }
        pts <- pts[-hull[which.min(crit)],]
    }
    out <- pts[chull(pts),]
    attr(out, "centre") <- polycentre(out)
    attr(out, "area") <- polyarea(out)
    out
}

## Area and centre of a polygon

`polyarea` <- function(x)
{
    x <- rbind(x, x[1,])
    n <- nrow(x)
    if (n < 4) # two points or less: area 0
        return(0)
    else
        abs(sum(x[-n,1]*x[-1,2] - x[-1,1]*x[-n,2]))/2
}

`polycentre` <- function(x)
{
    x <- rbind(x, x[1,])
    n <- nrow(x)
    if (n < 4)
        return(colMeans(x[-n,, drop = FALSE]))
    xy <- x[-n,1]*x[-1,2] - x[-1,1]*x[-n,2]
    A <- sum(xy)/2
    xc <- sum((x[-n,1] + x[-1,1]) * xy)/A/6
    yc <- sum((x[-n,2] + x[-1,2]) * xy)/A/6
    structure(c(xc, yc), names = colnames(x))
}

### Similar to peelhull but for ellipses: remove point with largest
### Mahalanobis distance, update and remove next point, and finally
### return the enclosing ellipse

#' @param pts Coordinates of points, a two-column matrix
#' @param keep Proportion of points kept
#'
#' @importFrom stats cov mahalanobis predict
#' @importFrom cluster ellipsoidhull volume
#' @rdname peelhull
#' @export
`peelellipse` <-
    function(pts, keep = 0.9, criterion = c("area", "mahalanobis"))
{
    stopifnot(ncol(pts) == 2, keep <= 1, keep > 0)
    criterion <- match.arg(criterion)
    ndrop <- as.integer(nrow(pts) * (1 - keep))
    for (k in seq_len(ndrop)) {
        del <-
            switch(criterion,
                   "mahalanobis" = which.max(mahalanobis(pts, colMeans(pts),
                                                         cov(pts))),
                   "area" = which.min(sapply(seq_len(nrow(pts)), function(i)
                                             volume(ellipsoidhull(pts[-i,])))))
            pts <- pts[-del,, drop=FALSE]
    }
    ell <- ellipsoidhull(pts)
    out <- predict(ell)
    attr(out, "centre") <- ell$loc
    attr(out, "area") <- volume(ell)
    out
}
