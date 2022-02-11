#' @importFrom stats rbeta
#' @importFrom vegan designdist
#' @export
`bayesjaccard` <-
    function(x, expected = FALSE)
{
    if (expected)
        designdist(x, "(b+c+1)/(a+b+c+2)", terms = "binary", abcd=TRUE)
    else
        designdist(x, "rbeta(length(a), b+c+1, a+1)", terms = "binary",
                   abcd = TRUE)
}

#' @importFrom stats hclust rbeta
#' @importFrom vegan designdist ordilabel
#' @export
`clsupport` <-
    function (x, n=1000, method="average", softmatch = FALSE, plot = TRUE, ...)
{
    h0 <- hclust(bayesjaccard(x, expected = TRUE), method)
    m0 <- clsets(h0)
    supp <- s <- numeric(nrow(h0$merge)-1)
    for(i in seq_len(n)) {
        d <- bayesjaccard(x)
        h <- hclust(d, method)
        m <- clsets(h)
        if (softmatch) {
            ## Find the Jaccard similarities of each m0 & m sets and
            ## take the highest as the value of "soft match". aabc
            ## (a+b + a+c) is the sum of species numbers of two sets,
            ## and abc the size of union of sets (number of units in
            ## pooled sets, a+b+c).
            aabc <- outer(sapply(m0, length), sapply(m, length), "+")
            abc <- matrix(0, length(m), length(m))
            for(i in seq_along(m))
                abc[,i] <- sapply(m0, function(x) length(unique(c(m[[i]], x))))
            abc <- aabc/abc - 1 # Jaccard similarity = a/abc
            s <- apply(abc, 1, max)
        } else {
        for (i in seq_along(supp))
            s[i] <- any(sapply(m, function(x) all(identical(m0[[i]], x))))
        }
        supp <- supp + s
    }
    if (plot) {
        ## with softmax supp is a float
        if (softmatch) {
            supp <- round(1000 * supp/n)
        }
        plot(h0, ...)
        ordilabel(h0, "internal", labels=supp)
    }
    supp
}

#' @export
`clsets` <-
    function(hclus)
{
    m <- hclus$merge
    memb <- list()
    sets <- list()
    for (j in seq_len(nrow(m) - 1)) {
        for (k in 1:2)
            sets[[k]] <- if (m[j, k] < 0)
                             -m[j,k]
                         else
                             memb[[m[j,k]]]
        memb[[j]] <- sort(c(sets[[1]], sets[[2]]))
    }
    memb
}

#############
### NMDS  ###
#############

#' @importFrom stats fitted
#' @importFrom vegan metaMDS monoMDS procrustes
#' @export
`bjNMDS` <-
    function(x, n = 100, trymax = 500, maxit = 1000, smin = 1e-4,
             sfgrmin = 1e-7, sratmax = 0.999999, parallel = 2, trace=FALSE,
             ...)
{
    ## Expected ordination
    d0 <- bayesjaccard(x, expected=T)
    m0 <- metaMDS(d0, trymax = trymax, maxit = maxit, smin = smin,
                  sfgrmin = sfgrmin, sratmax = sratmax, parallel = parallel,
                  trace = trace, ...)
    ## random Jaccard ordinations
    rscore <- array(dim = c(dim(m0$points), n))
    for (i in seq_len(n)) {
        d <- bayesjaccard(x)
        m <- monoMDS(d, m0$points, maxit = maxit, smin = smin,
                     sfgrmin = sfgrmin, sratmax = sratmax, ...)
        rscore[,,i] <- fitted(procrustes(m0$points, m$points), truemean=FALSE)
    }
    m0$rscores <- rscore
    class(m0) <- c("bjnmds", class(m0))
    m0
}

#' @importFrom graphics polygon
#' @importFrom grDevices col2rgb rgb
#' @export
`bjpolygon` <-
    function(x, keep = 0.9, kind = c("hull", "ellipse"), col="gray",
             alpha = 127, observed = TRUE, ...)
{
    if (!inherits(x, "bjnmds"))
        stop("needs bayesjaccard ordination object")
    kind <- match.arg(kind)
    xarr <- x$rscores
    x0 <- x$points
    dims <- dim(xarr)
    nobs <- dims[1]
    nsam <- dims[3]
    ## handle colours
    if (alpha < 1)
        alpha <- round(255 * alpha)
    if (is.factor(col))
        col <- as.numeric(col)
    cols <- rgb(t(col2rgb(col)), alpha = alpha, maxColorValue = 255)
    cols <- rep(cols, length = nobs)
    ## draw polygons
    for (i in seq_len(nobs)) {
        poly <- switch(
            kind,
            "hull" = peelhull(t(xarr[i,,]), keep = keep,
                              criterion = "distance"),
            "ellipse" = peelellipse(t(xarr[i,,]), keep = keep)
        )
        if (kind == "hull" && observed)
            poly <- peelhull(rbind(poly, x0[i,]), keep = 1)
        polygon(poly, col = cols[i], border = NA)
    }
}

#' @importFrom graphics points segments
#' @importFrom vegan ordilabel
#' @export
`bjstars` <-
    function(x, keep = 0.9, col="gray", type = c("t", "p", "n"), ...)
{
    if (!inherits(x, "bjnmds"))
        stop("needs bayesjaccard ordination object")
    xarr <- x$rscores
    x0 <- x$points
    nobs <- nrow(x0)
    nsam <- dim(xarr)[3]
    nkept <- ceiling(nsam * keep)
    if (nkept == nsam)
        kept <- rep(TRUE, nsam)
    col <- rep(col, length = nobs)
    for (i in seq_len(nobs)) {
        if (nkept < nsam) {
            dist <- colSums((xarr[i,,] - x0[i,])^2)
            kept <- rank(dist) <= nkept
        }
        segments(x0[i,1], x0[i,2], xarr[i,1,kept], xarr[i,2,kept], col = col[i],  ...)
    }
    switch(type,
           "n" = invisible(),
           "p" = points(x0, ...),
           "t" = ordilabel(x0, ...)
           )
    invisible()
}
