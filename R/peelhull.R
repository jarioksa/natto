#'
#' Convex Hull Enclosing a Given Proportion of Points
#'
#' Function finds a small convex hull enclosing a given proportion of
#' points. The function works by removing points from the hull one by
#' one. The points are selected either so that the area of the
#' remaining hull is reduced as much as possible at every step (de
#' Smith, Goodrich & Longley 2007), or removing the point that has the
#' maximal total distance to all remaining points. The function only
#' works in 2D.
#'
#' @encoding UTF-8
#'
#' @param pts Coordinates of points, a two-column matrix
#' @param keep Proportion of points kept
#' @param criterion Criterion to remove a point on the
#'     hull. \code{"area"} removes the point that reduces the area of
#'     the new hull as much as possible, and \code{"distance"} removes
#'     the point that has the largest total distance to all points on
#'     and within the hull.
#'
#' @details This is preliminary work (but without guarantee of
#'     progress). However, the current function is such that it can be
#'     easily plugged into \code{\link[vegan]{ordihull}}
#'     (\CRANpkg{vegan} package).
#'
#' The algorithms are naïve, and stepwise removal of single points
#' does not guarantee smallest possible final hull. Although the area
#' reduction algorithm was found in literature (de Smith et al. 2007),
#' the distance criterion often gives smaller final convex hulls.
#'
#' @author Jari Oksanen
#'
#' @references de Smith, M.J., Goodchild, M.F. & Longley,
#'     P.A. (2007). \emph{Geospatial analysis: A comprehensive guide
#'     to principles, techniques and software tools}. Matador.
#'
#' @return A two-column matrix of coordinates of peeled hull defining
#'     a closed convex hull that can be plotted with
#'     \code{\link[graphics]{polygon}}.
#'
#' @importFrom grDevices chull
#'
#' @export
`peelhull` <-
    function(pts, keep = 0.9, criterion = c("area", "distance"))
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
               "area" = polyarea(pts[chull(pts[-hull[i],]),]),
               "distance" = -sum(sweep(pts, 2, pts[hull[i],])^2))
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
    function(pts, keep = 0.9)
{
    stopifnot(ncol(pts) == 2, keep <= 1, keep > 0)
    ndrop <- as.integer(nrow(pts) * (1 - keep))
    for (k in seq_len(ndrop)) {
        del <- which.max(mahalanobis(pts, colMeans(pts), cov(pts)))
        pts <- pts[-del,, drop=FALSE]
    }
    ell <- ellipsoidhull(pts)
    out <- predict(ell)
    attr(out, "centre") <- ell$loc
    attr(out, "area") <- volume(ell)
    out
}
