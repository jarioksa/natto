#' @param x Fitted result object.
#' @param \dots Other parameters to functions.
#' @importFrom graphics points
#' @rdname humpfit
#' @export
`points.humpfit` <-
    function(x, ...)
{
    points(x$x, x$y, ...)
    invisible()
}
