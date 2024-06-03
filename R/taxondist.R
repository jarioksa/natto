### Clarke's taxonomic distance was suggested for vegan in
### https://github.com/vegandevs/vegan/issues/430. This function is
### based on a function of a response to that request.

#' Clarke's Taxonomic Community Dissimilarity
#'
#' Taxonomic dissimilarity takes into account related species when
#' there is no exact species match in compared communities. This
#' function returns generalizations of Bray-Curtis and Kulczynski
#' indices (Clarke _et al._ 2006).
#'
#' @details When a species occurs only in one of two compared
#'     communities, the species will increase dissimilarity of two
#'     communities with one species. However, if the other communities
#'     contains related species, the increase will be less depending
#'     on the relatedness of species (Clarke _et al._ 2006). The
#'     degree of relatedness is defined by taxonomic or other taxon
#'     dissimilarities.
#'
#'     The current function generalizes Bray-Curtis and Kulczynski
#'     indices (see \code{\link{canneddist}}) for binary
#'     (presence-absence) data by taking the dissimilarity to most
#'     closely related species as the increase of dissimilarity,
#'     instead of 1 of basic index. If species dissimilarities have
#'     any values over 1, they will be scaled by observed
#'     maximum. Alternatively, user can specify maximum species
#'     dissimilarity (argument `dmax`) for scaling, and if it is
#'     lower than data maximum, higher values will be truncated to
#'     scaled value 1. This allows imposing stricter concept of
#'     relatedness so that, say, taxa in different orders will be
#'     regarded as completely unrelated.
#'
#'     Although the function was proposed and is named as taxonomic
#'     dissimilarity, the function can be used with phylogenetic,
#'     functional or trait dissimilarities.
#'
#'     Function \code{\link{distrao}} provides an alternative index
#'     that can handle quantitative data and produce indices related
#'     to Euclidean distance.
#'
#' @seealso \code{\link{distrao}} is an alternative taxonomic distance
#'     measure. The index is closely related to Clarke's taxonomic
#'     diversity indices which are available in \pkg{vegan} function
#'     \code{\link[vegan]{taxondive}}. \pkg{Vegan} function
#'     \code{\link[vegan]{taxa2dist}} can be used to find taxonomic
#'     dissimilarities from classification table.
#' @examples
#'
#' if (require(vegan)) {
#' data(dune, dune.phylodis)
#' ## phylogenetic data, but regard lineaages completely distict
#' ## beyond K/T (K/Pg) age limit
#' d <- taxondist(dune, dune.phylodis, dmax = 65.2)
#' polarord(d)
#' }
#'
#' @param x Community data; will be treated as binary presence/absence
#'     matrix.
#' @param d Taxonomic, phylogenetic or other dissimilarities among
#'     species (columns of \code{x}).
#' @param dmax Scale dissimilarities by \code{dmax} (if \code{dmax >
#'     max(d)}) or truncate dissimilarities at \code{dmax} (if
#'     \code{dmax < max(d)}).
#' @param method Type of returned dissimilarity index. Gamma is
#'     generalized Bray-Curtis, and Theta generalized Kulczynski
#'     index.
#'
#' @return Clarke's taxonomic dissimilarity index as defined in
#'     \code{method}.
#'
#' @references Clarke, K.R., Somerfield, P.J. & Chapman,
#'     M.G. (2006). On resemblance measures for ecological studies,
#'     including taxonomic dissimilarities and a zero-adjusted
#'     Bray-Curtis coefficient for denuded
#'     assemblages. _J. Exp. Marine Biol. & Ecol._ 330, 55-80.
#'
#' @importFrom stats as.dist
#' @export
`taxondist` <-
    function (x, d, method = c("gamma", "theta"), dmax)
{
    method <- match.arg(method)
    x <- as.matrix(x)
    x <- ifelse(x > 0, 1, 0)
    d <- as.dist(d)
    ## scale or truncate by dmax
    if (!missing(dmax)) {
        d <- d/dmax
        if (max(d) > 1)
            d[d > 1] <- 1
    } else if (max(d) > 1) {
        d <- d/max(d)
    }
    d <- as.matrix(d)
    if (NCOL(d) != NCOL(x))
        stop("Number of columns do not match in 'x' and 'd'")
    N <- NROW(x)
    dis <- matrix(0, N, N)
    for(j in 1:(N-1)) {
        for(i in (j+1):N) {
            crosstd <- (outer(x[i,], x[j,]) * d)[x[i,] > 0, x[j,] > 0,
                                                       drop = FALSE]
            min1 <- apply(crosstd, 1, min)
            min2 <- apply(crosstd, 2, min)
            dis[i,j] <- switch(method,
                               "gamma" = mean(c(min1, min2)),
                               "theta" = (mean(min1) + mean(min2))/2)
        }
    }
    dis <- as.dist(dis)
    attr(dis, "call") <- match.call()
    attr(dis, "method") <- paste("clarke", method, sep=".")
    attr(dis, "Labels") <- rownames(x)
    attr(dis, "maxdist") <- 1
    dis
}
