#' @param xlab,ylab Axis labels.
#' @param lwd Line width.
#' @param l.col,p.col Line and point colours.
#' @param type Type of ype of \code{plot}: \code{"p"} for observed
#'     points, \code{"l"} for fitted lines, \code{"b"} for both, and
#'     \code{"n"} for only setting axes.
#' @importFrom graphics plot points lines
#' @rdname humpfit
#' @export
`plot.humpfit` <-
    function(x, xlab="Biomass", ylab="Species Richness", lwd=2, l.col="blue",
             p.col = 1, type="b", ...)
{
    plot(x$x, x$y, xlab = xlab, ylab = ylab, type="n", ...)
    if (type == "b" || type == "p")
        points(x, col = p.col, ...)
    if (type == "b" || type == "l")
        lines(x, lwd = lwd, col = l.col, ...)
    invisible()
}
