`diverclust` <-
    function (x, renyi = 1, equalize = TRUE, beta = TRUE, hill = FALSE,
          trace = TRUE, ...) 
{
    require(vegan) || stop("I need vegan")
    x <- as.matrix(x) ## huge speed-up over data frame
    ## equalize: divide each row by a scaling factor that allows
    ## getting arithmetic averages for baseline of beta diversity and
    ## pool rows; the scaling factor depends on 'renyi'
    if (equalize) {
        if (renyi==1)
           x <- decostand(x, "total")
        else if(renyi==2)
           x <- decostand(x, "norm", MARGIN=1)
        else if (renyi==Inf)
           x <- decostand(x, "max", MARGIN=1)
        else if (renyi==0)
           x <- decostand(x, "pa")
        else {
           tmp <- (rowSums(x^renyi))^(1/renyi)
           x <- sweep(x, 1, tmp, "/")
        }
    }
    cli <- seq_len(nrow(x))
    id <- -cli
    cnt <- rep(1, nrow(x))
    merge <- matrix(NA, nrow=nrow(x)-1, ncol=2)
    height <- rep(NA, nrow(x)-1)
    dh <- matrix(NA, nrow(x), nrow(x))
    ## for beta diversity we need 'base' of alpha diversities
    if (beta)
       base <- renyi(x, scale = renyi, hill=hill)
    else
       base <- numeric(nrow(x))
    ## Inial pooled (beta) diversities of all n*(n-1)/2 pairs of sites.
    if (trace)
        cat("Getting pooled diversity of all pairs of sites: this may take time\n")
    for(i in 2:length(cli))
        for(j in 1:(i-1)) {
            dh[i,j] <- renyi(x[i,]+x[j,], scales = renyi, hill = hill)
            dh[i,j] <- dh[i,j] - (base[i] + base[j])/2
    }
    ## Cluster
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
    if(hill) {
        dm <- switch(as.character(renyi),
                     "0" = "Species Richness",
                     "1" = "exp Shannon",
                     "2" = "Simpson",
                     "Inf" = "Berger-Parker",
                     paste("Hill", as.character(renyi)))
    } else {
        dm <- switch(as.character(renyi),
                     "0" = "log Species Richness",
                     "1" = "Shannon",
                     "2" = "log Simpson",
                     "Inf" = "log Berger-Parker",
                     paste("Renyi", as.character(renyi)))
    }
    if(beta)
        dm <- paste(dm, "beta")

    out <- list(merge = merge, 
                height = height,
                order=hclustMergeOrder(merge),
                labels = rownames(x),
                dist.method = dm,
                call = match.call(),
                method = "diverclust")
    class(out) <- c("hclust", "diverclust")
    out
}

### Internal vegan function to get the 'order' from a merge matrix of
### an hclust tree

`hclustMergeOrder` <-
    function(merge)
{
    ## Get order of leaves with recursive search from the root
    order <- numeric(nrow(merge)+1)
    ind <- 0
    ## "<<-" updates data only within hclustMergeOrder, but outside
    ## the visit() function.
    visit <- function(i, j) {
        if (merge[i,j] < 0) {
            ind <<- ind+1
            order[ind] <<- -merge[i,j]
        } else {
            visit(merge[i,j], 1)
            visit(merge[i,j], 2)
        }
    }
    visit(nrow(merge), 1)
    visit(nrow(merge), 2)
    return(order)
}
