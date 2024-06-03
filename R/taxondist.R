### Clarke's taxonomic distance was suggested for vegan in
### https://github.com/vegandevs/vegan/issues/430. This function is
### based on a function of a response to that request.

#' Return Gammaplus and Thetaplus Elements of Clarke's Taxonomic
#' Community Dissimilarity
#'
#' @param x Community data; will be treated as binary presence/absence
#'     matrix.
#' @param d Taxonomic, phylogenetic or other dissimilarities among
#'     species (columns of \code{x}).
#' @param dmax Scale dissimilarities by \code{dmax} (if \code{dmax >
#'     max(d)}) or truncate dissimilarities at \code{dmax} (if
#'     \code{dmax < max(d)}).
#' @param method Type of returned dissimilarity index.
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
