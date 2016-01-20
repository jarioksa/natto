"iap" <-
    function (comm, freq.min=20, permutations=1000) 
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
        q <- mean(spno[comm[,k]>0])
        out[i,2] <- q 
        out[i,1] <- freq[k]
        for (j in 1:permutations)
            sim[j] <- mean(sample(spno, freq[k], prob=spno))
        out[i,3] <- mean(sim)
        out[i,4:5] <- quantile(sim, c(0.05, 0.95))
        p <- if (q < median(sim)) sum(q < sim) else sum(q > sim)
        out[i,6] <- min(1, 2*(1 - p/permutations))
    }
    class(out) <- "iap"
    out
}

"plot.iap" <-
    function (x, ...) 
{
    plot.default(x, ...)
    i <- order(x[,1])
    matlines(x[i,1], x[i,3:5], lty=c(1,2,2), ...)
}

"print.iap" <-
    function (x, ...) 
{
    printCoefmat(x, ...)
    invisible(x)
}

"summary.iap" <-
    function (object, ...) 
{
    x <- object[object[,6] <= 0.1,]
    i <- order(x[,2])
    x <- x[i,]
    class(x) <- "iap"
    x
}

