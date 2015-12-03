#' Information Clusterings of Ecological Communities
#'
#' Function performs a cluster analysis based on information analysis
#' as defined by Lambert et al. (1966). The current implementation is
#' based on Legendre & Legendre (2012).
#'
#' @param x Community data
#' @param delta Use increase in information (Delta I) instead of the
#' value of information shown in the tree when merging clusters as
#' recommended by Lambert et al. (1966, see Legendre & Legendre 2012
#' for details). If \code{TRUE}, some \code{\link{hclust}} may not
#' work correctly with the result.
#'
#' @references Williams, W.T., Lambert, J.M. & Lance,
#' G.N. (1966). Multivariate methods in plant ecology. V. Similarity
#' analyses and information-analysis. \emph{J. Ecol.} 54, 427--445.
#'
#' Legendre, P. & Legendre, L. (2012). \emph{Numerical Ecology.} 3rd
#' English Ed., Elsevier.
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
    dis[lower.tri(dis)] <- -2 * log(0.5) *
        designdist(x, "A+B-2*J", "binary", name = "information")
    adj <- numeric(N)
    merge <- matrix(NA, nrow = N - 1, ncol = 2)
    height <- rep(NA, N - 1)
    id <- -seq_len(N)
    w <- rep(1, N)
    for (lev in 1:(N-1)) {
        ## pick minimum distance
        if (delta)
            g <- which.min(dis - outer(adj, adj, "+"))
        else
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
    class(out) <- c("hclust", "infoclust")
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
