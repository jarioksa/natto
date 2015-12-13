#' beta and gamma Diversity Distance
#'
#' The distance between two sampling units is the increase in
#' diversity when these are merged together, or the beta diversity
#' between two SUs. Alternatively, the function can find the total
#' diversity or distance from the origin which is known as gamma
#' diversity for two SUs. Any Renyi diversity can be used.
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
#' @return A dissimilarity object inheriting from \code{\link{dist}}.
#' @author Jari Oksanen
#'
#' @importFrom vegan decostand renyi
#' @importFrom stats as.dist
#' 
#' @export
`diverdist` <-
    function (x, renyi = 1, equalize = TRUE, beta = TRUE, hill = FALSE) 
{
    x <- as.matrix(x) ## huge speed-up over data frame
    ## equalize
    if (equalize) {
        x <- renyiEqualize(x, renyi = renyi)
    }
    dh <- matrix(NA, nrow(x), nrow(x))
    ## for beta diversity we need 'base' of alpha diversities
    if (beta)
        alpha <- renyi(x, scale = renyi, hill=hill)
    else
        alpha <- numeric(nrow(x))
    ## pooled (beta) diversities of all n*(n-1)/2 pairs of sites.
    dh[lower.tri(dh)] <-
        do.call(c, lapply(seq_len(nrow(x)-1), function(i)
            renyi(sweep(x[-seq_len(i),,drop =FALSE], 2, x[i,], "+"),
                  scale = renyi, hill = hill)))
    if (beta)
        dh <- dh - outer(alpha, alpha, "+")/2
    ## make to 'dist'
    dh <- as.dist(dh)
    ## method name
    if(hill) {
        dm <- switch(as.character(renyi),
                     "0" = "Species Richness",
                     "1" = "exp Shannon",
                     "2" = "Simpson",
                     "Inf" = "Berger-Parker",
                     paste("Hill", as.character(renyi)))
    } else {
        dm <- switch(as.character(renyi),
                     "0" = "log Species Richness",
                     "1" = "Shannon",
                     "2" = "log Simpson",
                     "Inf" = "log Berger-Parker",
                     paste("Renyi", as.character(renyi)))
    }
    if(beta)
        dm <- paste(dm, "beta")
    attr(dh, "method") <- dm
    attr(dh, "alpha") <- alpha
    attr(dh, "Labels") <- rownames(x)
    attr(dh, "call") <- match.call()
    dh
}

### internal function for equalizing community data: divide each row
### by a scaling factor that allows getting arithmetic averages for
### baseline of beta diversity and pool rows; the scaling factor
### depends on 'renyi'

`renyiEqualize` <-
    function(x, renyi = 1)
{
    if (renyi==1)
        decostand(x, "total")
    else if(renyi==2)
        decostand(x, "norm", MARGIN=1)
    else if (renyi==Inf)
        decostand(x, "max", MARGIN=1)
    else if (renyi==0)
        decostand(x, "pa")
    else {
        tmp <- (rowSums(x^renyi))^(1/renyi)
        sweep(x, 1, tmp, "/")
    }
}
