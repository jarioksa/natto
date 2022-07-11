#' Orthogonalization via QR Decomposition
#'
#' Function orthogonalizes vector \code{x} against all columns of
#' matrix \code{x0}.
#'
#' @param x0 Matrix against which orthogonalization is calculated.
#' @param x vector to be orthogonalized.
#'
#' @export
`ortho`  <-
    function(x0, x)
{
    qr.resid(qr(x0), x)
}
