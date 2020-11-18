#' @importFrom stats cov2cor
#' @importFrom vegan bstick
#'
#' @export
`resassoc.cca` <-
    function(x, rank, ...)
{
    if (!inherits(x, "cca"))
        stop("input must be a constrained ordination object from vegan")
    if (is.null(x$CA$eig))
        stop("object 'x' does not have residual ordination")
    ## use bstick to get the rank if not given
    if (missing(rank)) {
        bs <- bstick(m$CA$rank, m$CA$tot.chi)
        rank <- min(which(bs > m$CA$eig)) - 1
    }
    ## Residual associations can be estimated only when rank > 1, else
    ## return identity matrix
    if (rank < 2) {
        r <- diag(NROW(m$CA$v))
        dimnames(r) <- list(rownames(m$CA$v), rownames(m$CA$v))
    } else {
        eig <- sqrt(m$CA$eig[seq_len(rank)])
        x <- m$CA$v[, seq_len(rank)] %*% diag(eig, nrow=rank)
        r <- cov2cor(tcrossprod(x))
    }
    attr(r, "rank") <- rank
    r
}
