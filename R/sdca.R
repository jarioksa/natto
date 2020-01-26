#' Smoothly Detrended Correspondence Analysis
#'
#' Function \code{sdca} is similar to \code{\link[vegan]{decorana}},
#' but instead of detrending by segments it uses \code{\link{loess}}
#' for smooth non-linear detrending.
#'
#' @param Y Input data.
#' @param pairwise Detrend axis \eqn{k} separately for each previous
#'     axis in order \eqn{1 \dots k-1 \dots 1}. This only influences
#'     axes 3 and 4.
#' @param \dots Other arguments (passed to \code{\link{loess}}).

#' @importFrom stats loess residuals runif weighted.mean
#' @importFrom vegan eigengrad wascores

#' @export
`sdca` <-
    function(Y, pairwise = FALSE, ...)
{
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
                if (pairwise && k > 2) {
                   for(kk in c(1:(k-1), (k-2):1))
                       u <- residuals(loess(u ~ U[,kk], weights = r, ...))
                }
                else
                    u <- residuals(lo <- loess(u ~ U[,1:(k-1)], weights = r,
                                           normalize = FALSE,  ...))
                plot(u ~ U[,1], cex=0.3)
                points(U[,1], fitted(lo), pch=16, col=2)
            }
            vprime <- wascores(u, Y)
            eig <- abs(sum(c*v*vprime))
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
    colnames(V) <- paste0("DCA", seq_len(ncol(V)))
    colnames(U) <- paste0("DCA", seq_len(ncol(U)))
    rownames(V) <- colnames(Y)
    rownames(U) <- rownames(Y)
    origin <- apply(U, 2, weighted.mean, r)
    eig0 <- eigengrad(V, t(Y))
    out <- list(evals.decorana = EIG, evals = eig0, rproj = U, cproj = V,
                origin = origin, call = match.call())
    class(out) <- c("sdca", "decorana")
    out
}

## print method

#' @export
`print.sdca` <-
    function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n")
    cat(deparse(x$call), "\n\n")
    cat("Detrended correspondence analysis with loess smoothers\n\n")
    print(rbind(Eigenvalues = x$evals, `Decorana values` = x$evals.decorana),
          digits = digits)
    cat("\n")
    invisible(x)
}
