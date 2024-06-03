### standardize data matrix 'x' with distances 'd' so Simpson
### diversity will be Rao's quadratic entroy and Euclidean distances
### will be (transformed) Rao dissimilarity: x <- x %*% (1-d)^1/2
#' Standardize Data to Yield Rao Satistics
#'
#' Function \code{raostand} modified data so that it is compatible to
#' other functions, and Rao's quadratic entropy and Rao distances can
#' be directly found from the standarized data.
#'
#' @details
#'
#' Function \code{raostand} standardizes data similarlty as implicitly
#' done in \code{qrao} and \code{raodist} when \code{propx =
#' TRUE}. For standardized data \code{Z}, quadratic entropy is found
#' as \code{1 - rowSums(Z^2)}, and Rao distances can be found via
#' Euclidean distances of \code{Z}. The standardized data allows
#' calculating any generic community dissimilarity, and using \code{Z}
#' in \code{\link[vegan]{rda}} allows performing phylogenitically
#' constrained RDA. The standardization does not preserve absences,
#' but zero abundances will be boosted to positive values when the
#' sampling unit has related species.  More details can be found in
#' vignette.
#'
#' @examples
#' ## Rao standardization
#' ## Phylogenetic diversity
#' data(dune, dune.phylodis, dune.env)
#' Z <- raostand(dune, dune.phylodis)
#' all.equal(1 - rowSums(Z^2), qrao(dune, dune.phylodis),
#'      check.attributes = FALSE)
#' ## Phylogenetic distance
#' all.equal(dist(Z), distrao(dune, dune.phylodis, method="euclidean"),
#'      check.attributes = FALSE)
#' ## plot with standardized values
#' tabasco(Z, hclust(dist(Z)), hclust(dune.phylodis))
#' ## phylogenetic polar ordination
#' pol <- polarord(vegdist(Z))
#' sppscores(pol) <- Z
#' plot(pol)
#' ## Phylogenetically constrainted RDA
#' mod <- rda(raostand(dune, dune.phylodis, propx = FALSE) ~ Management + Moisture,
#'     dune.env)
#' anova(mod, by = "margin")
#' plot(mod, scaling = "sites")
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
#' @rdname qrao
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
    if (anyNA(d))
        stop("missing values are not accepted")
    if (!missing(dmax)) {
        d <- d/dmax
        if (max(d) > 1)
            d[d > 1] <- 1
    }
    else if (max(d) > 1)
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
    if (any(abs(x) < TOL))
        x[abs(x) < TOL] <- 0
    attr(x, "dimnames") <- dn
    x
}
