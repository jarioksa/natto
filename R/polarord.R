### Polar Ordination A.K.A. Bray-Curtis Ordination according to Bruce
### McCune & Jim Grace (2002), Analysis of Ecological Communities,
### Chapter 17

`polarord` <-
    function(d, k=2)
{
    ## Create a zero matrix of ordination scores
    axes <- matrix(0, attr(d, "Size"), k)
    colnames(axes) <- paste0("PO", seq_len(k))
    rownames(axes) <- attr(d, "Labels")
    ## Iterate through dimensions 1..k
    for(dim in seq_len(k)) {
        m <- as.matrix(d)
        ## Find the first endpoint (variance method)
        p1 <- which.max(apply(m, 2, var))
        ## Find the second endpoint (regression method)
        p2 <- which.min(cov(m[p1,],m))
        ## Project all points between these endpoints
        axes[,dim] <- (m[p1,p2]^2  + m[p1,]^2 - m[p2,]^2)/2/m[p1,p2]
        ## Find the residual dissimilarities, with a guard against
        ## negative residuals in semimetric indices
        d <- sqrt(pmax(d^2 - dist(axes[,dim])^2,0))
    }
    axes
}
