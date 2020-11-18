#' Residual Species Associations in Constrained Ordination
#'
#' Function finds residual species association after fitting the
#' constraints (and conditions) in constrained ordination
#' (\code{\link[vegan]{cca}}, \code{\link[vegan]{rda}}).
#'
#' @param x Constrained ordination result from
#'     \code{\link[vegan]{cca}} or \code{\link[vegan]{rda}}.
#' @param rank Number of unconstrained ordination axes that are used
#'     to compute the species associations. If missing, the rank is
#'     chosen by brokenstick distribution
#'     (\code{\link[vegan]{bstick}}).
#'
#' @return Correlation-like association matrix between species with
#'     added argument \code{"rank"} of the unconstrained data used to
#'     compute associations.
#'
#' @examples
#' library(vegan)
#' data(dune, dune.env)
#' mod <- cca(dune ~ A1 + Moisture, dune.env)
#' resassoc.cca(mod)
#'
#' @importFrom stats cov2cor
#' @importFrom vegan bstick
#'
#' @export
`resassoc.cca` <-
    function(x, rank)
{
    if (!inherits(x, "cca"))
        stop("input must be a constrained ordination object from vegan")
    if (inherits(x, c("dbrda","capscale")))
        stop("distance-based methods do not have species scores")
    if (is.null(x$CA$eig))
        stop("object 'x' does not have residual ordination")
    ## use bstick to get the rank if not given
    if (missing(rank)) {
        bs <- bstick(x$CA$rank, x$CA$tot.chi)
        rank <- min(which(bs > x$CA$eig)) - 1
    }
    ## Residual associations can be estimated only when rank > 1, else
    ## return identity matrix
    if (rank < 2) {
        r <- diag(NROW(x$CA$v))
        dimnames(r) <- list(rownames(x$CA$v), rownames(x$CA$v))
    } else {
        eig <- sqrt(x$CA$eig[seq_len(rank)])
        x <- x$CA$v[, seq_len(rank)] %*% diag(eig, nrow=rank)
        r <- cov2cor(tcrossprod(x))
    }
    attr(r, "rank") <- rank
    r
}
