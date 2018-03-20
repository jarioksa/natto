#' Convex Hull Enclosing the Given Proportion of Points
#'
#' Function finds a small convex hull enclosin a given proportion of
#' points. The function works by removing point from the hull that
#' gives the largest reduction in the area of remaining hull until a
#' desired number of points are removed (de Smith, Goodrich & Longley
#' 2007). The function only works in 2D.
#'
#' @param pts Coordinates of points, a two-column matrix
#' @param keep Proportion of points kept
#'
#' @details This is preliminary work (but without guarantee of
#'     progress). However, the current function is such that it can be
#'     easily plugged into \code{\link[vegan]{ordihull}}
#'     (\CRANpkg{vegan} package).
#'
#' @author Jari Oksanen
#'
#' @references de Smith, M.J., Goodchild, M.F. & Longley,
#'     P.A. (2007). \emph{Geospatial analysis: A comprehensive guide
#'     to principles, techniques and software tools}. Matador.
#'
#' @return A two-column matrix of coordinates of peeled hull defining
#'     a closed convex hull (last and first lines are the same).
#'
#' @importFrom grDevices chull
#'
#' @export
`peelhull` <-
    function(pts, keep = 0.9)
{
    stopifnot(ncol(pts) == 2, keep <= 1, keep > 0)
    ndrop <- as.integer(nrow(pts) * (1-keep))
    if (ndrop == 0)
        return(chull(pts))
    ## remove one point and select a new hull with smallest area
    for (k in seq_len(ndrop)) {
        hull <- chull(pts)
        area <- numeric(length(hull))
        for (i in seq_along(hull)) {
            area[i] <- polyarea(pts[chull(pts[-hull[i],]),])
        }
        pts <- pts[-hull[which.min(area)],]
    }
    pts[chull(pts),]
}

#" Find area of a polygon

`polyarea` <- function(x) {
    x <- rbind(x, x[1,])
    n <- nrow(x)
    x
    if (n < 4)
        return(0)
    else
        abs(sum(x[-n,1]*x[-1,2] - x[-1,1]*x[-n,2]))/2
}
