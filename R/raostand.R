### standardize data matrix 'x' with distances 'd' so Simpson
### diversity will be Rao's quadratic entroy and Euclidean distances
### will be (transformed) Rao dissimilarity: x <- x %*% (1-d)^1/2
#' @export
`raostand`<-
    function(x, d)
{
    TOL <- sqrt(.Machine$double.eps)
    x <- as.matrix(x)
    dn <- attr(x, "dimnames")
    ## distances d to similarity matrix
    if (max(d) > 1)
        d <- d/max(d)
    d <- as.matrix(d)
    d <- 1 - d
    ## check that diagonals are 1
    if (any(abs(diag(d) - 1) > TOL))
        stop("'d' is not a valid dissimilarity object")
    ## eigen decomposition
    e <- eigen(d)
    ## all eigenvalues should be positive
    if (any(e$values < -TOL))
        stop("dissimilarities 'd' do not define Euclidean transformation")
    k <- e$values > TOL
    ## matrix squareroot
    vec <- e$vectors[, k, drop=FALSE]
    ev <- e$values[k]
    d <- vec %*% (sqrt(ev) * t(vec))
    ## transform
    x <- x %*% d
    attr(x, "dimnames") <- dn
    x
}
