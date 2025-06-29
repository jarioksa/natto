#' Rao's Quadratic Entropy and Dissimilarity
#'
#' Rao's quadratic entropy (Rao 1982) is a generalization of Simpson
#' (or Gini-Simpson) diversity index to a situation where species are
#' non-independent. The species dependences are expressed as
#' dissimilarities, where 1 means independent species, and 0 a
#' completely aliased species. Rao's distance (1982) is a similar
#' generalization for distances between sampling units with
#' non-independent species. Typically species distances are based on
#' phylogenetics, taxonomy or species traits, and these measures are
#' called phylogenetic or functional diversities and distances.
#'
#' @details Rao's quadratic entropy is \eqn{H_i = \sum_{j} \sum_{k}
#' p_{ij} p_{ik} d_{jk}}{Q[i] = sum(k) sum(j) p[ij] * p[ik] * d[jk]},
#' where \eqn{i} is the index of the sampling unit, \eqn{j} and
#' \eqn{k} are indices of species, \eqn{p} is the proportion of
#' species in the sampling unit, and \eqn{d} are distances among
#' species. The distances should be scaled to range \eqn{0...1}, and
#' they are divided by the observed maximum if this exceeds
#' 1. Alternatively, the distances are divided by argument \code{dmax}
#' instead of data maximum. Distances that are shorter than
#' \code{dmax} are truncated to the maximum value. The square roots of
#' distances should be Euclidean, but this is not verified. They are
#' Euclidean if there are no negative eigenvalues in the principal
#' coordinates analysis, and \pkg{ade4} package has function
#' \code{\link[ade4]{is.euclid}} for a canned test. If all distances
#' are 1, species are independent and \code{qrao} will return Simpson
#' diversity.
#'
#' Function \code{distrao} finds distances based on quadratic entropy
#' for sampling units. Rao (1982) suggested distance \eqn{d_{ij} =
#' H_{ij} - \frac{1}{2}(H_i + H_j)}{d[ij] = H[ij] - 0.5*(H[i]+H[j])},
#' where \eqn{H_i}{H[i]} and \eqn{H_j}{H[j]} are Rao entropies for
#' sites \eqn{i} and \eqn{j} and \eqn{H_{ij}}{H[ij]} is the similar
#' entropy evaluated so that the species proportions \eqn{p} are from
#' different sampling units. Rao (1982) called this as Jensen
#' distance, and it is half of the squared Euclidean distance. The
#' Euclidean distance can also be requested. In addition, Rao (1982)
#' suggested a standardized distance that is based on logarithms of
#' elements \eqn{H}.
#'
#' @param x Community data. The diversity is estimated for the rows of
#'     the matrix.
#' @param d Dissimilarities among species (columns). If some some
#'     dissimilarities are \eqn{>1}, these are divided with the
#'     maximum dissimilarity. If all \eqn{d = 1} or \code{d} is
#'     missing, \code{qrao} will return Simpson's index, and
#'     \code{qraodist} a basic dissimilarity index.
#' @param na.rm Should missing values be removed?
#' @param dmax Scale dissimilarities by \code{dmax} (if
#'     \code{dmax > max(d)}) or truncate dissimilarities at \code{dmax}
#'     (if \code{dmax < max(d)}).
#' @return \code{qrao} returns a vector of Rao's quadratic entropy
#'     values and \code{distrao} distances of class \code{"dist"}.
#'
#' @references Rao, C.R. (1982) Diversity and dissimilarity
#' coefficients: a unified approach. \emph{Theoretical Population
#' Biology} 21, 24--43.
#'
#' @seealso There are other implementations of this function in
#'     \R. Most notably functions \code{\link[ade4]{divc}} and
#'     \code{\link[ade4]{disc}} in \pkg{ade4}. However, these may
#'     square input dissimilarities and divide results by 2 depending
#'     on options.  Function \code{\link{taxondist}} provides Clarke's
#'     taxonomic dissimilarity.
#'
#' @examples
#' if (require(vegan)) {
#' data(dune, dune.phylodis)
#' qrao(dune, dune.phylodis)
#' tabasco(dune, hclust(distrao(dune, dune.phylodis)), hclust(dune.phylodis))
#' ## regard lineages completely distinct beyond K/T (K/Pg) limit
#' qrao(dune, dune.phylodis, dmax = 65.2)
#' }
#'
#' @importFrom vegan decostand
#'
#' @export
`qrao` <-
    function(x, d, na.rm = FALSE, dmax)
{
    ## handle community matrix
    x <- as.matrix(x)
    x <- decostand(x, "tot")
    ## dissimilarities among species columns. These should be in 0..1,
    ## where 1 imeans completely different species, and 0 completely
    ## identical.
    if (missing(d)) {
        d <- 1
    } else {
        d <- as.dist(d)
        ## scale d by dmax > max(d) or truncate d at dmax < max(d)
        if (!missing(dmax)) {
            d <- d/dmax
            if (max(d, na.rm = na.rm) > 1)
                d[d > 1] <- 1
        } else if (max(d, na.rm = na.rm) > 1)
            d <- d/max(d, na.rm = na.rm)
    }
    ## diversities
    n <- nrow(x)
    p <- ncol(x)
    ltri <- lower.tri(matrix(0, p, p))
    div <- numeric(n)
    for(i in seq_len(n))
        div[i] <- 2*sum(outer(x[i,], x[i,])[ltri] * d, na.rm = na.rm)
    div
}
#' @param method Distance measure to be used (see Details).
#' @importFrom stats as.dist
#' @rdname qrao
#' @export
`distrao` <-
    function(x, d, method = c("jensen", "euclidean", "standardized"), dmax)
{
    method <- match.arg(method)
    x <- as.matrix(x)
    if (missing(d)) {
        d <- diag(nrow = ncol(x))
    } else {
        d <- as.dist(d)
        ## scale or truncate by dmax (NB, no NA allowed)
        if (!missing(dmax)) {
            d <- d/dmax
            if (max(d) > 1)
                d[d > 1] <- 1
        } else if (max(d) > 1)
            d <- d/max(d)
        d <- as.matrix(d)
    }
    ## qrao found only diagonal elements, but now we need off-diagonal, too.
    H <- x %*% d %*% t(x)
    ## Distances
    if (method == "standardized")
        H <- log(H)
    out <- outer(diag(H), diag(H), "+")/2 - H
    if (method == "euclidean")
        out <- sqrt(2*out)
    out <- as.dist(out)
    attr(out, "method") <- paste("rao", method, sep = ".")
    attr(out, "call") <- match.call()
    attr(out, "maxdist") <- NA
    out
}
