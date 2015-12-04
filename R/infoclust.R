#' Information Clusterings of Ecological Communities
#'
#' Function performs hierarchical clustering of binary ecological
#' communities based on information analysis as defined by Williams et
#' al. (1966) and Lance & Williams (1966).
#'
#' @param x Community data
#' @param delta Use increase in information (\eqn{\Delta I}) instead of the
#' value of information (\eqn{I}) when merging clusters as
#' recommended by Williams et al. (1966) and Lance & Williams (1966).
#'
#' @details Function performs information analysis of binary
#' ecological communities (Williams et al. 1966, Lance & Williams
#' 1966). The current implementation is based on Legendre & Legendre
#' (2012).
#'
#' The information \eqn{I} of a collection of \eqn{N} sampling units
#' with \eqn{S} species is defined as
#'   \deqn{I = S N \log N - \sum_i^S a_i \log a_i + (N-a_i) \log (N-a_i)}{I = S*N*log(N) - sum(a*log(a) + (N-a)*log(N-a))}
#' where \eqn{a} is the frequency count of each species in the
#' collection. The method works by merging either the units that give
#' the lowest increase (\eqn{\Delta I} when \code{delta = TRUE}), or
#' the units that are most homogeneous (lowest \eqn{I} when
#' \code{delta = FALSE}). After merging sampling units or clusters,
#' the community data matrix is updated by actually merging the data
#' units and re-evaluating their information distance to all other
#' units. The information content of all non-merged clusters is \eqn{I
#' = 0}, and for clusters of several sampling units the constant
#' species (completely absent or always present) do not contribute to
#' the information. The largest increase in information is made by
#' species with 0.5 relative frequency, so that the analysis tries to
#' build clusters where species is either always present or always
#' absent. This often gives easily interpretable clusters.
#'
#' @references Williams, W.T., Lambert, J.M. & Lance,
#' G.N. (1966). Multivariate methods in plant ecology. V. Similarity
#' analyses and information-analysis. \emph{J. Ecol.} 54, 427--445.
#'
#' Lance, G.N. & Williams, W.T. (1966). Computer programs for
#' hierarchical polythetic classification (\dQuote{similarity
#' analyses}). \emph{Comp. J.} 9, 60--64.
#'
#' Legendre, P. & Legendre, L. (2012). \emph{Numerical Ecology.} 3rd
#' English Ed., Elsevier.
#'
#' @return Function returns an object of class \code{"infoclust"} that
#' inherits from \code{\link{hclust}}. It uses all \code{"hclust"}
#' methods, but some may fail or work in unexpected ways because the
#' analysis is not based on dissimilarities but on binary data matrix.
#'
#' @examples
#' ## example used to demonstrate the calculation of
#' ## information analysis by Legendre & Legendre (2012, p. 372).
#' data(pond)
#' cl <- infoclust(pond)
#' plot(cl, hang = -1)
#' ## Lance & Williams suggest a limit below which clustering is
#' ## insignificant and should not be interpreted
#' abline(h=qchisq(0.95, ncol(pond)), col=2)
#'
#' @importFrom vegan designdist
#' 

#' @export
`infoclust` <-
    function(x, delta = TRUE)
{
    requireNamespace("vegan") || stop("requires vegan package")
    x <- ifelse(x > 0, 1, 0)
    N <- nrow(x)
    ## shortcut to pairwise 'infodist'
    dis <- matrix(NA, N, N)
    dis[lower.tri(dis)] <- 2 * log(2) *
        designdist(x, "A+B-2*J", "binary", name = "information")
    adj <- numeric(N)
    merge <- matrix(NA, nrow = N - 1, ncol = 2)
    height <- rep(NA, N - 1)
    id <- -seq_len(N)
    w <- rep(1, N)
    for (lev in 1:(N-1)) {
        ## pick minimum distance
        if (delta) # use increase Delta I
            g <- which.min(dis - outer(adj, adj, "+"))
        else       # use the level of heterogeneity I
            g <- which.min(dis)
        g1 <- row(dis)[g]
        g2 <- col(dis)[g]
        ## update tree
        merge[lev,] <- c(id[g1], id[g2])
        height[lev] <- dis[g1,g2]
        id[g1] <- lev
        ## update distances
        if (delta)
            adj[g1] <- adj[g2] <- dis[g1,g2]
        x[g1,] <- x[g1,] + x[g2,]
        w[g1] <- w[g1] + w[g2]
        dis[g2,] <- NA
        dis[,g2] <- NA
        for (j in 1:N) {
            if (is.na(dis[max(g1,j), min(g1,j)]) || g1 == j)
                next
            d <- infodist(x[g1,] + x[j,], w[g1] + w[j])
            dis[max(g1,j), min(g1,j)] <- d
        }
    }
    ## height and merge should be ordered in hclust objects, but with
    ## 'delta' the clusters are not formed at the height order. We
    ## need to reorder these so that other R functions know how to
    ## handle the result.
    if (delta && is.unsorted(height)) {
        o <- order(height)
        oo <- order(o)
        for (i in 1:nrow(merge))
            for (j in 1:2)
                if (merge[i,j] > 0) # merge a (reordered?) cluster
                    merge[i,j] <- oo[merge[i,j]]
        merge <- merge[o,]
        height <- height[o]
    }
    out <- list(merge = merge,
                height = height,
                order = vegan:::hclustMergeOrder(merge),
                labels = rownames(x),
                dist.method = "information",
                call = match.call(),
                method = "infoclust")
    class(out) <- c("infoclust", "hclust")
    out
}

## unexported internal function
`infodist` <-
    function(a, n)
{
    length(a) * n * log(n) -
        sum(a * log(a), na.rm=TRUE) -
        sum((n-a)*log(n-a), na.rm = TRUE)
}
