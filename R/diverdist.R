##' beta and gamma Diversity Distance
#'
#' The distance between two sampling units is the increase in
#' deversity when these are merged together, or the beta diversity
#' between fused units. \code{diverdist} finds such distances for all
#' pairs of sampling units using any Rényi diversity.
#'
#' @return \code{diverdist} returns a dissimilarity object inheriting
#'     from \code{\link{dist}}.
#'
#' @importFrom vegan decostand renyi
#' @importFrom stats as.dist
#'
#' @examples
#' ## increase in Shannon diversity when pooling two sites
#' d <- diverdist(spurn, renyi=1, equalize="no")  # raw data
#' de <- diverdist(spurn, renyi=1, equalize="total")
#' plot(d, de, xlab = "Shannon distance, no equalization",
#'             ylab = "Shannon distance, equalized")
#'
#' @rdname diverclust
#'
#' @importFrom vegan decostand
#' @export
`diverdist` <-
    function (x, renyi = 1, equalize = c("renyi", "total", "no"),
              beta = TRUE, hill = FALSE)
{
    x <- as.matrix(x) ## huge speed-up over data frame
    ## equalize
    if (is.logical(equalize))
        equalize <- if (equalize) "renyi" else "no"
    equalize <- match.arg(equalize)
    if (equalize == "total") {
        x <- decostand(x, "total")
    } else if (equalize == "renyi") {
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

#' @importFrom vegan decostand

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
