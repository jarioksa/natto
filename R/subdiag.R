#' Extract Subdiagonal from Dissimilarities or from a Square Matrix
#'
#' @param x Square matrix or a \code{\link{dist}} object.
#' 
#' @export
`subdiag` <-
    function(x)
{
    if (inherits(x, "data.frame"))
        x <- as.matrix(x)
    if (is.matrix(x)) {
        if (ncol(x) != nrow(x))
            stop("only square matrix accepted")
        k <- row(x) == col(x)+1
    }
    else if(inherits(x, "dist")) {
        N <- attr(x, "Size")
        if (N <= 2)
            k <- 1
        else
            k <- c(1, cumsum((N-1) : 2) + 1)
    } else {
        stop("invalid input")
    }
    x[k]
}

