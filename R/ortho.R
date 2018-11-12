#' Gram-Schmidt Orthogonalization
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
    for(i in 1:ncol(x0))
        x <- x - drop(crossprod(x0[,i],x)/crossprod(x0[,i],x0[,i])) * x0[,i]
    x
}
