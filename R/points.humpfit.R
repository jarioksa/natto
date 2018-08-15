#' @importFrom graphics points
#' @export
`points.humpfit` <-
    function(x, ...)
{
    points(x$x, x$y, ...)
    invisible()
}
