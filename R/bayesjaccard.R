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
    function (x, n=1000, method="average", plot = TRUE, ...)
{
    h0 <- hclust(bayesjaccard(x, expected = TRUE), method)
    m0 <- clsets(h0)
    supp <- s <- numeric(nrow(h0$merge)-1)
    for(i in seq_len(n)) {
        d <- bayesjaccard(x)
        h <- hclust(d, method)
        m <- clsets(h)
        for (i in seq_along(supp))
            s[i] <- any(sapply(m, function(x) all(identical(m0[[i]], x))))
        supp <- supp + s
    }
    if (plot) {
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

