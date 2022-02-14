#' @importFrom stats rbeta
#' @importFrom vegan designdist
#' @export
`bayesjaccard` <-
    function(x, method = c("rbeta", "jaccard1"))
{
    method <- match.arg(method)
    switch(method,
           "jaccard1" =
               designdist(x, "(b+c+1)/(a+b+c+2)", terms = "binary", abcd=TRUE),
           "rbeta" =
               designdist(x, "rbeta(length(a), b+c+1, a+1)", terms = "binary",
                          abcd = TRUE)
           )
}

#' @importFrom stats hclust rbeta
#' @importFrom vegan designdist ordilabel
#' @export
`clsupport` <-
    function (x, n=1000, method="average", softmatch = FALSE, plot = TRUE, ...)
{
    h0 <- hclust(bayesjaccard(x, method="jaccard1"), method)
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
    d0 <- bayesjaccard(x, method="jaccard1")
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
#' @importFrom grDevices adjustcolor col2rgb rgb
#' @export
`bjpolygon` <-
    function(x, keep = 0.9, kind = c("hull", "ellipse"), linetopoint = TRUE,
             col="gray", alpha = 127, observed = TRUE, ...)
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
        if (linetopoint) {
            cnt <- attr(poly, "centre")
            ## line is non-transparent
            segments(x0[i,1], x0[i,2], cnt[1], cnt[2],
                     col = adjustcolor(cols[i], 1))
        }
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


#############
### dbRDA ###
#############

### Somewhat trickier to implement than NMDS. (1) dbRDA is an
### eigenvector method, and if we rotate, eigenvalues would change,
### and we may need to skip rotation, but fix the axis reflection.
### (2) The "expected" model has one negative eigenvalue, but rbeta
### models can have several and variable numbers of negative
### eigenvalues. (3) LC scores can be rotation invariant when the real
### ranks are equal, but can become kinkier with different real ranks
### (but we perhaps we should not rotate). (4) Probably we should not
### allow any adjustment against negative eigenvalues as these are
### data-set dependent and destroy the beauty of the distance (excpet
### sqrt.dis?).

#' @importFrom vegan dbrda eigenvals
#' @export
`bjdbrda` <-
    function(formula, data, n=100, ...)
{
    environment(formula) <- environment()
    m0 <- dbrda(formula = formula, data = data, distance="jaccard1",
                dfun = bayesjaccard, ...)
    if (is.null(m0$CCA))
        stop("only implememented for constrained analysis")
    naxes <- length(m0$CCA$eig)
    ## sample, but first collect only eigenvalues
    tot.chi <- numeric(n)
    ev <- matrix(NA, nrow = n, ncol = naxes)
    u0 <- m0$CCA$u
    u <- array(dim = c(dim(u0), n))
    wa <- array(dim = c(dim(u0), n))
    bp <- array(dim = c(dim(m0$CCA$biplot), n))
    cn <- array(dim = c(dim(m0$CCA$centroids), n))
    for (i in 1:n) {
        m <- dbrda(formula, data, distance="rbeta", dfun = bayesjaccard, ...)
        tot.chi[i] <- m$tot.chi
        ev[i,] <- eigenvals(m, model = "constrained")
        ## no Procrustes rotation, but check reflected axes
        nreal <- seq_len(ncol(m$CCA$u))
        rev <- sign(colSums(u0[, nreal] * m$CCA$u))
        u[,nreal,i] <- m$CCA$u %*% diag(rev)
        wa[,nreal,i] <- m$CCA$wa %*% diag(rev)
        bp[,nreal,i] <- m$CCA$biplot %*% diag(rev)
        cn[,nreal,i] <- m$CCA$centroids %*% diag(rev)
    }
    BJ <- list("tot.chi" = tot.chi, "eig" = ev, "u" = u, "wa" = wa,
               "biplot" = bp, "centroids" = cn)
    m0$BayesJaccard <- BJ
    class(m0) <- c("bjdbrda", class(m0))
    m0
}
