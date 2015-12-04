`diverclust` <-
    function (x, renyi = 1, equalize = TRUE, beta = TRUE, hill = FALSE,
          trace = TRUE, ...) 
{
    x <- as.matrix(x) ## huge speed-up over data frame
    ## equalize here, but no more in diverdist
    if (equalize)
        x <- renyiEqualize(x, renyi = renyi)
    cli <- seq_len(nrow(x))
    id <- -cli
    cnt <- rep(1, nrow(x))
    merge <- matrix(NA, nrow=nrow(x)-1, ncol=2)
    height <- rep(NA, nrow(x)-1)
    dh <- matrix(NA, nrow(x), nrow(x))
    if (trace)
        cat("Getting pooled diversities for all n*(n-1)/2 pairs of sites...")
    dis <- diverdist(x, renyi = renyi, equalize = FALSE, beta = beta,
                     hill = hill)
    dh[lower.tri(dh)] <- dis
    ## base: alpha diversities (or 0 if beta = FALSE)
    base <- attr(dis, "alpha")
    ## Cluster
    if (trace)
        cat(" done\nStarting clustering\n")
    for(lev in 1:(nrow(x)-1)) {
        ## g1 and g2 are indices of smallest (possibly tied) pooled
        ## (beta) diversity
        g1 <- row(dh)[which.min(dh)]
        g2 <- col(dh)[which.min(dh)]
        ## basic hclust statistics: merge & height
        merge[lev,] <- c(id[g1], id[g2])
        height[lev] <- dh[g1, g2]
        ## weights for pooling data
        w <- c(cnt[g1], cnt[g2])/(cnt[g1]+cnt[g2])
        ## In the following we basically pool data ('x') to the row
        ## 'g1' and make 'g2' NA. First we need counts
        cnt[g1] <- cnt[g1] + cnt[g2]
        ## update 'id': for single points this is a negative index,
        ## and for clusters the number of level of last merge
        id[g1] <- lev
        ## pool data and average alpha diversity
        x[g1,] <- w[1]*x[g1,] + w[2]*x[g2,]
        dh[g2,] <- NA
        dh[,g2] <- NA
        cli[g2] <- NA
        base[g1] <- w[1] * base[g1] + w[2] * base[g2]
        if (trace)
            cat(lev,": ", merge[lev,], " at ", height[lev], "\n")
        ## update pooled diversities for recently pooled data
        for(i in 1:nrow(x)) {
           if(is.na(cli[i]) || i == g1) next
           w <- c(cnt[g1], cnt[i])/(cnt[g1] + cnt[i])
           tmp <- renyi(w[1]*x[g1,] + w[2]*x[i,], scales = renyi, hill=hill)
           dh[max(g1,i), min(g1,i)] <- tmp - w[1]*base[g1] - w[2]*base[i] 
        }
    }

    out <- list(merge = merge, 
                height = height,
                order = vegan:::hclustMergeOrder(merge),
                labels = rownames(x),
                dist.method = attr(dis, "method"),
                call = match.call(),
                method = "diverclust")
    class(out) <- c("hclust", "diverclust")
    out
}
