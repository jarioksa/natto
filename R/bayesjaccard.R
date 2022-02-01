clsupport <-
    function (x, n=1000, method="average", usesets = FALSE, plot = TRUE, ...)
{
    h0 <- hclust(designdist(x, "(b+c+1)/(a+b+c+2)", term="bin", abcd=TRUE),
                 method)
    m0 <- if (usesets) clsets(h0) else clmembers(h0)
    supp <- s <- numeric(nrow(h0$merge)-1)
    for(i in seq_len(n)) {
        d <- designdist(x, "rbeta(length(a), b+c+1, a+1)",
                        term="bin", abcd=TRUE)
        h <- hclust(d, method)
        if (usesets) {
            m <- clsets(h)
            for (i in seq_along(supp))
                s[i] <- any(sapply(seq_along(supp), function(j)
                    setequal(m0[[i]], m[[j]])))
        } else {
            m <- clmembers(h)
            s <- sapply(seq_len(ncol(m)), function(i)
                any(apply(sweep(m, 1, m0[,i], "=="), 2, all)))
        }
        supp <- supp + s
    }
    if (plot) {
        plot(h0, ...)
        ordilabel(h0, "internal", labels=supp)
    }
    supp
}


clmembers <-
    function (hclus)
{
    m <- hclus$merge
    memb <- matrix(FALSE, nrow(m)+1, nrow(m)-1)
    for (j in seq_len(nrow(m) - 1)) {
        for (k in 1:2) {
            if(m[j,k] < 0)
                memb[-m[j,k], j] <- TRUE
            else
                memb[,j] <- memb[,j] | memb[,m[j,k]]
        }
    }
    memb
}

clsets <-
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
        memb[[j]] <- do.call(union, sets)
    }
    memb
}

