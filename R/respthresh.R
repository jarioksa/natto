### originally suggested to Dave Roberts as `dichtomod` in email on
### 30/4/2020.

`respthresh` <-
    function(object)
{
    fv <- fitted(object)
    fvuniq <- sort(unique(fv))
    y <- object$y
    fam <- object$family
    n <- length(fvuniq)
    dev <- numeric(n)
    mu0 <- mean(y)
    for(i in seq_len(n)) {
        out <- fv >= fvuniq[i]
        mu <- tapply(y, out, mean)
        dev[i] <- sum(fam$dev.resids(mu, mu0, table(out)))
    }
    hit <- which.max(dev)
    out <- list(threshold = fvuniq[hit], bestdev = dev[hit],
                cutoff = fvuniq, deviance = dev)
    class(out) <- "respthresh"
    out
}
