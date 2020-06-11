### originally suggested to Dave Roberts as `dichtomod` in email on
### 30/4/2020.

`respthresh` <-
    function(object)
{
    fv <- fitted(object)
    fvuniq <- sort(unique(fv))
    y <- object$y
    fam <- object$family
    wt <- weights(mod)
    n <- length(fvuniq)
    dev <- numeric(n)
    mu0 <- mean(y)
    for(i in 1:n) {
        out <- fv >= fvuniq[i]
        mu <- tapply(y, out, mean)
        dev[i] <- sum(fam$dev.resids(mu, mu0, table(out)))
    }
    list(cutoff = fvuniq, deviance = dev)
}
