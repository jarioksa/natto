#' @inheritParams points.humpfit
#' @param segments Number of segments used for lines.
#' @importFrom stats fitted predict
#' @importFrom graphics lines
#' @rdname humpfit
#' @export
`lines.humpfit` <-
    function(x, segments=101,  ...)
{
    mass <- x$x
    if (!is.null(segments) && segments > 0) {
        mass <- seq(min(mass), max(mass), length=segments)
        fv <- predict(x, newdata = mass)
    }
    else {
        i <- order(mass)
        fv <- fitted(x)
        mass <- mass[i]
        fv <- fv[i]
    }
    lines(mass, fv, ...)
    invisible()
}
