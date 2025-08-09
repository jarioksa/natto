#' Smoothly Detrended Correspondence Analysis
#'
#' Function \code{sdca} is similar to \code{\link[vegan]{decorana}},
#' but instead of detrending by segments it uses \code{\link{loess}}
#' for smooth non-linear detrending.
#'
#' Detrended correspondence analysis (DCA) as implemented in
#' \code{\link[vegan]{decorana}} tries to remove all systematic biases
#' between ordination axes by detrending later axes against previous
#' ones (Hill & Gauch 1980). DCA uses an ingenious method of
#' detrending by axis segments to allow removing non-linear
#' dependencies. Detrending means taking residuals against smoothed
#' segment means as the new ordination scores for the current axis. It
#' has been suggested that abrupt changes at segment borders can
#' cause some problems in DCA. The current function replaces segmented
#' detrending with detrending against smooth \code{\link{loess}}
#' functions. However, in many cases this changes the results little
#' from the original detrending by segments.
#'
#' The detrending for axes 3 and 4 is perfomed either using all
#' previous axes simultaneously in \code{\link{loess}} (default) or if
#' \code{pairwise=TRUE} performing sequential separate detrending
#' against each previous axis separately so that both the first and
#' last detrending are performed against axis 1 similarly as in the
#' original \code{decorana}. For axis 3 the sequence is against axes
#' 1, 2, 1 and for axis 4 against axes 1, 2, 3, 2, 1.
#'
#' The \code{\link[vegan]{decorana}} software made several other
#' innovations than detrending. Rescaling of axes is often more
#' influential in application than the actual detrending, like you can
#' see by using \code{\link[vegan]{decorana}} without rescaling.
#'
#' @param Y Input data.
#' @param iweigh Downweight rare species.
#' @param pairwise Detrend axis \eqn{k} separately for each previous
#'     axis in order \eqn{1 \dots k-1 \dots 1}. This only influences
#'     axes 3 and 4.
#' @param monitor Turn on graphical monitoring of detrending for each
#'     axis.
#' @param \dots Other arguments (passed to \code{\link{loess}}).

#' @importFrom stats loess residuals runif weighted.mean
#' @importFrom vegan downweight eigengrad wascores
#'
#' @return Function returns a subset of \code{\link[vegan]{decorana}}
#'     result object and can use many \code{decorana} methods (such as
#'     \code{plot}).
#'
#' @examples
#'
#' data(spurn)
#' mod <- sdca(spurn)
#' plot(mod, display = "sites")
#' plot(mod, display = "species", optimize = TRUE)
#' ## compare against rdecorana without rescaling
#' mod0 <- rdecorana(spurn, iresc = 0)
#' mod0
#' if (require(vegan)) {
#' plot(procrustes(mod0, mod, choices=1:2))
#' }

#' @references
#'
#'   Hill, M.O. and Gauch, H.G. (1980). Detrended correspondence
#'   analysis: an improved ordination technique. \emph{Vegetatio}
#'   \strong{42}, 47--58.

#' @export
`sdca` <-
    function(Y, iweigh = FALSE, pairwise = FALSE, monitor = FALSE, ...)
{
    EPS <- sqrt(.Machine$double.eps)
    CYCLES <- 200
    if (iweigh)
        Y <- downweight(Y, fraction=5)
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
                if (monitor) {
                    plot(u ~ U[,1], cex=0.3)
                    points(U[,1], fitted(lo), pch=16, col=2)
                }
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
                origin = origin, iweigh = iweigh, call = match.call())
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
    cat("Detrended correspondence analysis with loess smoothers.\n")
    if (x$iweigh)
        cat("Downweigthing of rare species from fraction 1/5.\n")
    cat("\n")
    print(rbind(Eigenvalues = x$evals, `Decorana values` = x$evals.decorana),
          digits = digits)
    cat("\n")
    invisible(x)
}
