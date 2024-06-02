### Clarke's taxonomic distance was suggested for vegan in
### https://github.com/vegandevs/vegan/issues/430. This function is
### based on a function of a response to that request.

#' Return Gammaplus and Thetaplus Elements of Clarke's Taxonomic
#' Community Dissimilarity
#'
#' @param x Community data; will be treated as binary presence/absence
#'     matrix.
#' @param d Taxonomic, phylogenetic or other dissimilarity
#'     matrix among species (columns of x).
#' @param dmax Truncate dissimilarities to \code{dmax}.
#' @param method Type of returned dissimilarity measure.
#'
#' @return Clarke's taxonomic dissimilarity index as defined in
#'     \code{type}.
#'
#' @importFrom stats as.dist
#' @export
`taxondist` <-
    function (x, d, method = c("gamma", "theta"), dmax)
{
    method <- match.arg(method)
    x <- as.matrix(x)
    x <- ifelse(x > 0, 1, 0)
    if (!missing(dmax) && any(d > dmax))
        d[d > dmax] <- dmax
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
    attr(dis, "maxdist") <- NA
    dis
}
