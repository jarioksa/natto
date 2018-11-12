posvector <-
    function (x, scale = FALSE)
{
    x <- scale(x, scale = scale)
    u <- matrix(NA, nrow(x), ncol(x))
    colnames(u) <- paste0("Sp", seq_len(ncol(u)))
    for(i in 1:ncol(u)) {
        xx <- cov(x)
        xx <- xx/sqrt(diag(xx))
        crit <- rowSums(xx^2)
        if (max(crit, na.rm = TRUE) < sqrt(.Machine$double.eps))
            break
        take <- which.max(crit)
        print(names(take))
        u[,i] <- x[,take]
        colnames(u)[i] <- names(take)
        x <- apply(x, 2, function(z) ortho(u[,i, drop=FALSE], z))
    }
    u[, !is.na(colSums(u))]
}
