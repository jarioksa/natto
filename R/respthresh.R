### originally suggested to Dave Roberts as `dichtomod` in email on
### 30/4/2020.

#' Binary Threshold for Continuous Fitted Response
#'
#' Species distribution models (SDM) often use methods that give
#' continuous fitted responses for the probability of species
#' occurrence. Some users want to find a threshold that splits these
#' continuous responses to binary outcomes of presence and absence
#' (Allouche, Tsoar & Cadmon 2006). This function finds a threshold
#' that maximizes the explained deviance and gives as true
#' approximation of the original continuous response as possible.
#'
#' @param object Fitted \code{\link{glm}} or \code{\link[mgcv]{gam}}
#'     object. The function was developed for binary observations
#'     responses fitted with \code{\link[stats]{binomial}} error
#'     distribution, but it may process other types of models with
#'     unspecified and untested results: proceed at your own risk.
#'
#' @return
#'
#' The function returns an object of class \code{"respthresh"} with
#' following elements:
#' \itemize{
#'   \item \code{threshold}: optimal splitting threshold of fitted values.
#'   \item \code{bestdeviance}: the explained deviance at \code{threshold}.
#'   \item \code{cutoff, deviance}: sorted possible thresholds (i.e.,
#'     estimated fitted values) and associated explained deviances.
#' }
#'
#' @references
#'
#' Allouche, O., Tsoar, A. & Kadmon, R. (2006) Assessing the
#' accurraccy of species distribution models: prevelance, kappa and
#' the true skill statistic (TSS). \emph{Journal of Applied Ecology}
#' 43, 1223--1232.
#'
#' @importFrom stats fitted
#' @export
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
    fit <- fv >= fvuniq[hit]
    vals <- tapply(y, fit, mean)
    fit <- vals[fit+1]
    out <- list(threshold = fvuniq[hit], bestdev = dev[hit],
                cutoff = fvuniq, deviance = dev, values = vals,
                fitted = fit)
    class(out) <- "respthresh"
    out
}

### method functions

#' @export
`print.respthresh` <-
    function(x, ...)
{
    cat("Best threshold", x$threshold, "\n")
    cat("average fitted values", x$values, "\n")
    cat("explained deviance", x$bestdev, "\n")
}

#' @export
`plot.respthresh` <-
    function(x, ...)
{
    plot(x$cutoff, x$deviance, xlab = "Response Cutoff",
         ylab = "Explained Deviance", type = "l", ...)
    abline(v = x$threshold, col = 2, ...)
}

