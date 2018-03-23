#' Clustering Based on Arrhenius Species-Area Model
#'
#' Function forms hierarchic clustering so that each fusion minimizes
#' the exponent \eqn{z} of Arrhenius species-area model.
#' #'
#' @param x Community data for clustering.
#'
#' @details The Arrhenius species area model defines the observed
#'     total diversity \eqn{\gamma}{gamma} as a function of average
#'     \eqn{\alpha}{alpha} diversity, number of sampling units \eqn{N}
#'     as \eqn{\gamma = \alpha N^z}{gamma = alpha*N^z}. The exponent
#'     \eqn{z} is a parameter of the steepness of the species-area
#'     -curve. When combining two clusters,
#'     \eqn{z=(\log(\gamma)-\log(\alpha))/\log(N)}{z = (log(gamma)-log(alpha))/log(N)}.
#'     The function combines clusters that have the lowest value of \eqn{z}.
#'
#' @return Function returns an \code{\link{hclust}} object.
#' @author Jari Oksanen.
#' @seealso \code{\link{hclust}} for cluster analysis and its support
#' methods, \code{\link{diverclust}} for related diversity clustering.
#'
#' @importFrom vegan designdist
#'
#' @export
`Azclust` <-
    function (x)
{
    x <- as.matrix(x)
    x <- ifelse(x > 0, 1, 0)
    cli <- seq_len(nrow(x))
    id <- -cli
    cnt <- rep(1, nrow(x))
    merge <- matrix(NA, nrow=nrow(x)-1, ncol=2)
    height <- rep(NA, nrow(x)-1)
    alpha <- rowSums(x)
    z <- matrix(NA, nrow(x), nrow(x))
    z[lower.tri(z)] <- designdist(x, "(log(gamma) - log(alpha))/log(2)",
                                  alphagamma = TRUE)
    ## Cluster
    for(lev in 1:(nrow(x)-1)) {
        ## g1 and g2 are indices of smallest (possibly tied) pooled
        ## (beta) diversity
        g1 <- row(z)[which.min(z)]
        g2 <- col(z)[which.min(z)]
        ## basic hclust statistics: merge & height
        merge[lev,] <- c(id[g1], id[g2])
        height[lev] <- z[g1, g2]
        ## In the following we basically pool data ('x') to the row
        ## 'g1' and make 'g2' NA. First we need counts
        w <- c(cnt[g1], cnt[g2])/(cnt[g1]+cnt[g2])
        cnt[g1] <- cnt[g1] + cnt[g2]
        ## update 'id': for single points this is a negative index,
        ## and for clusters the number of level of last merge
        id[g1] <- lev
        ## pool data and average alpha diversity
        x[g1,] <- pmax(x[g1,], x[g2,])
        z[g2,] <- NA
        z[,g2] <- NA
        cli[g2] <- NA
        alpha[g1] <- w[1] * alpha[g1] + w[2] * alpha[g2]
        ## update z for recently pooled data
        for(i in 1:nrow(x)) {
           if(is.na(cli[i]) || i == g1) next
           w <- c(cnt[g1], cnt[i])/(cnt[g1] + cnt[i])
           base <- w[1]*alpha[g1] + w[2]*alpha[i]
           gamma <- sum(pmax(x[g1,], x[i,]))
           tmp <- (log(gamma) - log(base))/log(cnt[g1]+cnt[i])
           z[max(g1,i), min(g1,i)] <- tmp
        }
    }


    out <- list(merge = merge,
                height = height,
                order = vegan:::hclustMergeOrder(merge),
                labels = rownames(x),
                call = match.call(),
                method = "diverclust")
    class(out) <- c("hclust", "diverclust")
    out
}
