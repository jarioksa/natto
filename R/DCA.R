`DCA` <-
    function(Y, ...)
{
    require(vegan) || stop("requires wascores function in vegan")
    EPS <- sqrt(.Machine$double.eps)
    CYCLES <- 200
    Y <- Y/sum(Y)
    r <- rowSums(Y)
    c <- colSums(Y)
    V <- matrix(0, ncol(Y), 4)
    U <- matrix(0, nrow(Y), 4)
    EIG <- numeric(4)
    v <- runif(ncol(Y))
    eig <- eig.old <- 0
    for (k in 1:4) {
        cycles <- 0
        repeat {
            v <- v - weighted.mean(v, w=c)
            v <- v/sqrt(sum(c*v*v))
            u <- wascores(v, t(Y))
            if (k > 1) {
                u <- residuals(lo <- loess(u ~ U[,1:(k-1)], weights = r,  ...))
                plot(u ~ U[,1], cex=0.3)
                points(U[,1], fitted(lo), pch=16, col=2)
            }
            vprime <- wascores(u, Y)
            eig <- sum(c*v*vprime)
            if ((tol <- abs(eig - eig.old)) < EPS || CYCLES <
                (cycles <- cycles+1)) {
                break
            } else {
                v <- vprime
                eig.old <- eig
            }
        }
        EIG[k] <- eig
        V[,k] <- v
        U[,k] <- wascores(v, t(Y))
        v <- runif(length(v))
        if (tol > EPS)
            warning("residual ", formatC(tol, digits=4), " larger than tolerance in axis ", k)
    }
    eig0 <- eigengrad(V, t(Y))
    list(DCA.eig = EIG, eig = eig0, u = U, v = V)
}
