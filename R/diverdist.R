`diverdist` <-
    function (x, renyi = 1, equalize = TRUE, beta = TRUE, hill = FALSE) 
{
    x <- as.matrix(x) ## huge speed-up over data frame
    ## equalize: divide each row by a scaling factor that allows
    ## getting arithmetic averages for baseline of beta diversity and
    ## pool rows; the scaling factor depends on 'renyi'
    if (equalize) {
        if (renyi==1)
            x <- decostand(x, "total")
        else if(renyi==2)
            x <- decostand(x, "norm", MARGIN=1)
        else if (renyi==Inf)
            x <- decostand(x, "max", MARGIN=1)
        else if (renyi==0)
            x <- decostand(x, "pa")
        else {
            tmp <- (rowSums(x^renyi))^(1/renyi)
            x <- sweep(x, 1, tmp, "/")
        }
    }
    cli <- seq_len(nrow(x))
    dh <- matrix(NA, nrow(x), nrow(x))
    ## for beta diversity we need 'base' of alpha diversities
    if (beta)
        alpha <- renyi(x, scale = renyi, hill=hill)
    else
        alpha <- numeric(nrow(x))
    ## pooled (beta) diversities of all n*(n-1)/2 pairs of sites.
    for(i in 2:length(cli))
        for(j in 1:(i-1)) {
            dh[i,j] <- renyi(x[i,]+x[j,], scales = renyi, hill = hill)
            dh[i,j] <- dh[i,j] - (alpha[i] + alpha[j])/2
        }
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
