#' beta and gamma Diversity Clustering
#'
#' \code{diverclust} forms hierarchic clustering so that each fusion
#' minimizes beta diversity or gamma diversity.
#'
#' @param trace Trace calculations. Either logical or an integer:
#' \code{trace=2} also traces merges.
#'
#' @details
#'
#' \code{diverclust} forms clusters so that pooled diversity or
#' its change to the baseline alpha diversity is minimized at each
#' level of clustering. The change of pooled diversity with respect to
#' the baseline alpha diversity is called beta diversity, and the
#' overall diversity is called gamma diversity. For beta diversity,
#' the clustering implies an additive partitioning.
#'
#' Function \code{diverdist} finds beta diversities (or optionally
#' gamma diversities) for pairs of sampling units. These are used as
#' starting values in \code{diverclust}, but after merging sampling
#' units into clusters, beta diversities are calculated for clusters
#' and cannot be found directly from these dissimilarities.
#'
#' Both \code{diverdist} and \code{diverclust} work by pooling
#' together sampling units and averaging species abundances. The units
#' may be equalized before pooling together to remove the effects of
#' variable totals in sampling units. Rényi diversities are based on
#' species proportions and most popular equalization is to use species
#' proportions by sampling units instead of raw data (argument
#' \code{equalize="total"}). With alternative \code{"renyi"}
#' equalization depends on the\code{renyi} scale. With \code{renyi=2}
#' (Simpson type) sums of squared proportions or norms are made equal,
#' and with \code{renyi=Inf}, species maxima are made equal in
#' compared units, and with \code{renyi=1} (Shannon type),
#' \code{"total"}. Equalization has no effect with \code{scale=0}
#' (Species richness type).
#'
#' For more detailed description, see vignette \code{"diverclust"}.
#'
#' @param x Community data.
#' @param renyi The scale of Renyi diversity index as defined in
#' \code{\link[vegan]{renyi}}. Value \code{0} gives diversity based on
#' species richness, \code{1} diversity based on Shannon entropy,
#' \code{2} diversity based on Simpson index, and \code{Inf} diversity
#' based on Berger-Parker index. Other non-negative values are also
#' allowed, but they may not define a well known standard diversity
#' index.
#' @param equalize Equalize data rows so that they can be meaningfully
#' pooled and averaged. If this is \code{FALSE}, raw data will be
#' used, and beta diversities may be negative. The equalization
#' depends on the value of \code{renyi}, and each row of \code{x} is
#' divided by \code{(rowSums(x^renyi))^(1/renyi)}. In borderline cases
#' \code{renyi=0} the data are presence-absence tranformed (but can be
#' left untransformed), and with \code{renyi=Inf} the rows are divided
#' by row maxima (see \code{\link[vegan]{decostand}}).
#' @param beta Use beta diversities: the average alpha diversity of
#' cluster members is subtracted from the pooled diversity. If this is
#' \code{FALSE},the clustering is based on pooled diversity, also
#' known as gamma diversity.
#' @param hill Use Hill numbers instead of Renyi diversity (see
#' \code{\link[vegan]{renyi}}). For \code{renyi = 0} these are species
#' richess values instead of their logarithms, and for \code{renyi =
#' 1} they are exponents of Shannon diversity instead of Shannon
#' diversities. In general, a Hill number is an exponent of Renyi
#' diversity. The Hill numbers may not be strictly additive, and beta
#' diversities may be negative (except with \code{renyi = 0}).
#'
#' @return \code{diverclust} returns an \code{\link{hclust}} object.
#' @author Jari Oksanen.
#' @seealso \code{\link{hclust}} for cluster analysis and its support
#'     methods, \code{\link[vegan]{renyi}} for estimating Renyi
#'     diversities and Hill numbers, and \code{\link[vegan]{adipart}}
#'     for related additive partitioning of beta diversity.
#'
#' @examples
#' data(spurn)
#' (cl <- diverclust(spurn))
#' plot(cl) # beta diversity based on Shannon diversity
#' (cl <- diverclust(spurn, renyi=0, hill=TRUE))
#' plot(cl) # increase in species richness when pooling sites
#'
#' @importFrom vegan decostand
#' @export
`diverclust` <-
    function (x, renyi = 1, equalize = c("renyi", "total", "no"),
              beta = TRUE, hill = FALSE, trace = FALSE)
{
    x <- as.matrix(x) ## huge speed-up over data frame
    ## equalize here, but no more in diverdist
    if (is.logical(equalize)) # for compatibility with old code
        equalize <- if (equalize) "renyi" else "no"
    equalize <- match.arg(equalize)
    if (equalize == "total")
        x <- decostand(x, "total")
    else if (equalize == "renyi")
        x <- renyiEqualize(x, renyi = renyi)
    cli <- seq_len(nrow(x))
    id <- -cli
    cnt <- rep(1, nrow(x))
    merge <- matrix(NA, nrow=nrow(x)-1, ncol=2)
    height <- rep(NA, nrow(x)-1)
    dh <- matrix(NA, nrow(x), nrow(x))
    if (trace)
        cat("Getting pooled diversities for all n*(n-1)/2 pairs of sites...")
    dis <- diverdist(x, renyi = renyi, equalize = FALSE, beta = beta,
                     hill = hill)
    dh[lower.tri(dh)] <- dis
    ## base: alpha diversities (or 0 if beta = FALSE)
    base <- attr(dis, "alpha")
    ## Cluster
    if (trace)
        cat(" done\nStarting clustering\n")
    for(lev in 1:(nrow(x)-1)) {
        ## g1 and g2 are indices of smallest (possibly tied) pooled
        ## (beta) diversity
        g1 <- row(dh)[which.min(dh)]
        g2 <- col(dh)[which.min(dh)]
        ## basic hclust statistics: merge & height
        merge[lev,] <- c(id[g1], id[g2])
        height[lev] <- dh[g1, g2]
        ## weights for pooling data
        w <- c(cnt[g1], cnt[g2])/(cnt[g1]+cnt[g2])
        ## In the following we basically pool data ('x') to the row
        ## 'g1' and make 'g2' NA. First we need counts
        cnt[g1] <- cnt[g1] + cnt[g2]
        ## update 'id': for single points this is a negative index,
        ## and for clusters the number of level of last merge
        id[g1] <- lev
        ## pool data and average alpha diversity
        x[g1,] <- w[1]*x[g1,] + w[2]*x[g2,]
        dh[g2,] <- NA
        dh[,g2] <- NA
        cli[g2] <- NA
        base[g1] <- w[1] * base[g1] + w[2] * base[g2]
        if (trace > 1)
            cat(lev,": ", merge[lev,], " at ", height[lev], "\n")
        ## update pooled diversities for recently pooled data
        for(i in 1:nrow(x)) {
           if(is.na(cli[i]) || i == g1) next
           w <- c(cnt[g1], cnt[i])/(cnt[g1] + cnt[i])
           tmp <- renyi(w[1]*x[g1,] + w[2]*x[i,], scales = renyi, hill=hill)
           dh[max(g1,i), min(g1,i)] <- tmp - w[1]*base[g1] - w[2]*base[i]
        }
    }

    out <- list(merge = merge,
                height = height,
                order = vegan:::hclustMergeOrder(merge),
                labels = rownames(x),
                dist.method = attr(dis, "method"),
                call = match.call(),
                method = "diverclust")
    class(out) <- c("hclust", "diverclust")
    out
}
