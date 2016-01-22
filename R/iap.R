#' Number of Companion Species (IAP)
#'
#' Function finds the number of companion species or the average
#' species richness of other species for each species in a
#' community. This is the quality index \eqn{Q} of Index of
#' Atmospheric Purity (IAP) in lichen bioindication (LeBlanc & De
#' Sloover 1970). The non-randomness of low or high \eqn{Q} values is
#' found by randomization.
#'
#' @details Index of Atmospheric Purity (IAP) is used in bioindication
#' with epiphytic lichens and bryophytes (LeBlanc & De Sloover
#' 1970). It derives species indicator scores (\eqn{Q}) as the number
#' of other species in sampling units where each focal species is
#' present, and then finds the IAP values for each sampling unit as
#' scaled weighted sum of species indicator values.
#'
#' Function \code{iap} finds the \eqn{Q} values for all species in a
#' community data set. This is a general measure of indicator value
#' for species richness and can well be used outside lichen
#' bioindication. The \eqn{Q} value is the average species richness in
#' sampling units where the species is present, excluding the species
#' itself from the richness. For rare species, \eqn{Q} is based on
#' small sample size, and is therefore more variable than for common
#' species. The \code{iap} function assesses the non-randomness
#' (\sQuote{significance}) of \eqn{Q} by taking random samples of the
#' same size as the frequency (number of occurrence) of the focal
#' species and finding the average richness (without the focal
#' species) in these samples. Because species are more likely to be
#' present in species-rich sampling units than in species-poor, the
#' random sampling uses the observed species richness (without the
#' focal species) as weights in random sampling. Testing is two-sided
#' and the number of greater or less random values is multiplied with
#' two. The observed value of \eqn{Q} is included in the random sample
#' of species richness values both in assessing the \eqn{p}-value and
#' in estimating the quantiles.
#' 
#' @references LeBlanc, S.C. & De Sloover, J. (1970) Relation between
#' industrialization and the distribution and growth of epiphytic
#' lichens and mosses in Montreal. \emph{Can. J. Bot.} 48, 1485--1496.
#' 
#' @param comm The community data frame.
#' @param freq.min Minimum number of occurrences for analysed species.
#' @param permutations Number of permutations to assess the randomized
#' number of companion species.
#'
#' @importFrom stats median quantile
#' @rdname iap
#' @export
`iap` <-
    function (comm, freq.min = 5, permutations = 999)
{
    spno <- rowSums(comm > 0)
    freq <- colSums(comm > 0)
    take <- freq >= freq.min
    ntake <- sum(take)
    idx <- seq(ncol(comm))[take]
    out <- matrix(NA, nrow=ntake, ncol=6)
    rownames(out) <- colnames(comm)[take] 
    colnames(out) <- c("Freq", "Q", "E(Q)", "5%", "95%", "Pr(Q = E(Q))")
    sim <- numeric(permutations)
    for (i in 1:ntake) {
        k <- idx[i]
        othersp <- spno - (comm[,k] > 0)
        q <- mean(othersp[comm[,k]>0])
        out[i,2] <- q 
        out[i,1] <- freq[k]
        for (j in 1:permutations)
            sim[j] <- mean(sample(othersp, freq[k], prob=othersp))
        out[i,3] <- mean(sim)
        out[i,4:5] <- quantile(c(sim, q), c(0.05, 0.95))
        p <- if (q < median(sim)) sum(q >= sim) else sum(q <= sim)
        out[i,6] <- min(1, (2*p+1)/(permutations+1))
    }
    class(out) <- "iap"
    out
}

#' @importFrom graphics matlines
#' @param x \code{iap} result object.
#' @param \dots Other arguments to the function.
#' @rdname iap
#' @export
`plot.iap` <-
    function (x, ...) 
{
    plot.default(x, ...)
    i <- order(x[,1])
    matlines(x[i,1], x[i,3:5], lty=c(1,2,2), ...)
}

#' @importFrom stats printCoefmat
`print.iap` <-
    function (x, ...) 
{
    printCoefmat(x, ...)
    invisible(x)
}

#' @param object \code{iap} result object.
#' @rdname iap
#' @export
summary.iap` <-
    function (object, ...) 
{
    x <- object[object[,6] <= 0.1,]
    i <- order(x[,2])
    x <- x[i,]
    class(x) <- "iap"
    x
}

