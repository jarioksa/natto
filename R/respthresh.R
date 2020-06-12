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
#'   \item \code{coefficients}: optimal splitting threshold of fitted values.
#'   \item \code{expl.deviance}: explained deviance
#'   \item \code{deviance}: residual deviance at \code{threshold}.
#'   \item \code{cutoff, devprofile}: sorted possible thresholds (i.e.,
#'     estimated fitted values) and associated explained deviances.
#'   \item \code{values}: fitted values below and above threshold.
#'   \item \code{fitted}: fitted values for each observation.
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
    out <- list(coefficients = fvuniq[hit], expl.deviance = dev[hit],
                deviance = mod$null.deviance-dev[hit],
                null.deviance = mod$null.deviance,
                orig.deviance = deviance(mod), formula = formula(mod),
                cutoff = fvuniq, devprofile = dev, values = vals,
                fitted = fit)
    class(out) <- "respthresh"
    out
}

### method functions

#' @export
`print.respthresh` <-
    function(x, ...)
{
    cat("Binary threshold for model\n")
    print(formula(x))
    cat("\nthreshold", x$coefficients, "\n")
    cat("average fitted values", x$values, "\n")
    cat("explained deviance", x$expl.deviance, "\n")
}

#' @export
`plot.respthresh` <-
    function(x, ...)
{
    plot(x$cutoff, x$devprofile, xlab = "Response Cutoff",
         ylab = "Explained Deviance", type = "l", ...)
    abline(v = x$coefficients, col = 2, ...)
}

#' @export
`summary.respthresh` <-
    function(object, ...)
{
    tmp <- c(object$null.deviance, deviance(object), object$orig.deviance)
    devtable <- cbind(tmp, c(NA, -diff(tmp)))
    dimnames(devtable) <- list(c("Null", "Binary", "Model"),
                               c("Deviance", "Delta Dev"))
    explained <- object$expl.deviance
    ok <- object$expl.deviance - object$devprofile <= qchisq(0.95, 1)
    oktable <- data.frame("Threshold" = object$cutoff[ok],
                          "Delta Dev" = object$expl.deviance - object$devprofile[ok],
                          check.names = FALSE)
    out <- list(devtable = devtable, oktable = oktable,
                formula = formula(object))
    class(out) <- "summary.respthresh"
    out
}

#' @export
`print.summary.respthresh` <-
    function(x, ...)
{
    cat("Binary threshold for model\n")
    print(formula(x))
    cat("\nTable of Deviance\n")
    printCoefmat(x$devtable, na.print="")
    cat("\nGood threshold values (the best with Delta Dev 0)\n")
    print(x$oktable)
}

