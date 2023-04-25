### standardize data matrix 'x' with distances 'd' so Simpson
### diversity will be Rao's quadratic entroy and Euclidean distances
### will be (transformed) Rao dissimilarity: x <- x %*% (1-d)^1/2
#' Standardize Data to Yield Rao Satistics
#'
#' Function standardizes data so that it gives Rao phylogenetic
#' diversity and Euclidean distance give Rao phylogenetic distance.
#'
#' @importFrom vegan decostand
#'
#' @param x Community data.
#' @param d Phylogenetic distances or other dissimilarities among
#'     species.
#' @param propx Standardize row of \code{x} to unit total. This will
#'     give standardization that is compatible with \code{qrao} and
#'     \code{distrao}.
#' @param dmax Scale dissimilarities by \code{dmax} (if \code{dmax >
#'     max(d)}) or truncate dissimilarities at \code{dmax} (if
#'     \code{dmax < max(d)}).
#'
#' @export
`raostand`<-
    function(x, d, propx = TRUE, dmax)
{
    TOL <- sqrt(.Machine$double.eps)
    x <- as.matrix(x)
    if (propx)
        x <- decostand(x, "tot")
    dn <- attr(x, "dimnames")
    ## distances d to similarity matrix
    d <- as.dist(d)
    if (!missing(dmax)) {
        d <- d/dmax
        if (max(d, na.rm = na.rm) > 1)
            d[d > 1] <- 1
    }
    else if (max(d, na.rm = na.rm) > 1)
        d <- d/max(d, na.rm = na.rm)
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
