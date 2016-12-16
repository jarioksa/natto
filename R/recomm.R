#' Rebuild Communities from Rarefied Pieces
#'
#' Community is divided into given number of pieces, and randomly
#' rarefied community is derived for each of these pieces using
#' \code{\link[vegan]{rrarefy}}. The final community is a sum of
#' rarefied pieces. The purpose is to mimick sampling variability of
#' the observed community. However, the process loses rare species,
#' and generated communities have lower species richness and diversity
#' than the original one.
#'
#' @return Randomized community data.

#' @param x communities to be rebuilt
#' @param pieces number of pieces used to rebuild the community
#'
#' @importFrom vegan rrarefy
#'
#' @author Jari Oksanen

#' @export
`recomm` <-
    function(x, pieces=4)
{
    x <- as.matrix(x)
    tot <- rowSums(x)
    splits <- round(tot/pieces)
    split1 <- tot - (pieces-1)*splits
    out <- rrarefy(x, split1)
    for(i in seq_len(pieces-1))
        out <- out + rrarefy(x, splits)
    out
}

#' Generate Communities as Random Poisson from Swan Expectations
#'
#' Poisson random numbers lose species as they introduce new zeros for
#' low expected values. They are also unable to compensate this by
#' generating new "unseen" species. To compensate this, we first use
#' Swan transformation with one pass to generate above-zero expected
#' values for unseen species. After generating the Swan probabilities
#' for missing species we standardize all species to their original
#' totals; this is a similar idea as the species coverage in the
#' Good-Turing model (Good 1953). Finally, a community is generated as
#' a Poisson random variate of adjusted observed abundances and Swan
#' probabilities. The process maintains average species richness and
#' diversity, but Poisson distribution probably underestimates the
#' real sampling variance.
#'
#' @author Jari Oksanen
#'
#' @return Randomized community matrix.
#'
#' @references Good, I.J. (1953) The population frequencies of species
#'     and the estimation of population parameters. Biometrika 40,
#'     237--264.

#' @param x community data of integer counts
#'
#' @importFrom vegan swan
#' @importFrom stats rpois
#'
#' @export
`swanrandom` <-
    function(x)
{
    x <- as.matrix(x)
    dims <- dim(x)
    dnam <- dimnames(x)
    coltot <- colSums(x)
    x <- swan(x, maxit = 1)
    x <- sweep(x, 2, coltot/colSums(x), "*")
    x <- rpois(prod(dims), x)
    dim(x) <- dims
    dimnames(x) <- dnam
    x
}
