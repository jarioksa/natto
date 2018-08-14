#' @importFrom graphics points
`points.humpfit` <-
    function(x, ...)
{
    points(x$x, x$y, ...)
    invisible()
}
