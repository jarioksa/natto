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
#' @details
#'
#' The function evaluates the deviance of the original continuous
#' model by setting a threshold at each fitted value and evaluating
#' the deviance explained with this threshold. The threshold is
#' inclusive: values at or above the threshold are regarded as being
#' hits. The deviance profile is usually jagged, and the method is
#' sensitive to single data points, and can give disjunct regions of
#' nearly equally good threshold points. Function \code{plot} will
#' display the deviance profile, and \code{summary} lists all
#' threshold that are nearly as good as the best point using
#' Chi-square distribution with 1 degree of freedom at \eqn{p=0.95} as
#' the criterion.
#'
#' The binary model has two values or average responses below and at
#' or above the threshold, but these values are not usually 0 and
#' 1. However, they are the values that maximize the explained
#' deviance of a two-value model.
#'
#' The \code{summary} also gives a Deviance table showing the original
#' Null deviance, the original model deviance, and between these the
#' residual deviance with the binary threshold, and the differences of
#' these deviances in the second model. The result object can also be
#' accessed with \code{\link{coef}} that returns the threshold,
#' \code{\link{deviance}} that returns the residual deviance, and
#' \code{\link{fitted}} that returns the fitted two values for each
#' original observation.
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
#'   \item \code{formula, null.deviance, orig.deviance}:
#'      \code{\link{formula}}, null deviance and deviance of the input
#'      model.
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
    if (object$family != "binomial" || length(unique(object$y)) > 2)
        warning("developed for binomial binary models: proceed at your own risk")
    fv <- fitted(object)
    fvuniq <- sort(unique(fv))
    y <- object$y
    fam <- object$family
    n <- length(fvuniq)
    dev <- numeric(n)
    mu0 <- mean(y)
    for(i in seq_len(n)) {
        above <- fv >= fvuniq[i]
        mu <- tapply(y, above, mean)
        dev[i] <- sum(fam$dev.resids(mu, mu0, table(above)))
    }
    hit <- which.max(dev)
    fit <- fv >= fvuniq[hit]
    vals <- tapply(y, fit, mean)
    fit <- vals[fit+1]
    out <- list(coefficients = fvuniq[hit], expl.deviance = dev[hit],
                deviance = object$null.deviance-dev[hit],
                null.deviance = object$null.deviance,
                orig.deviance = deviance(object), formula = formula(object),
                cutoff = fvuniq, devprofile = dev, values = vals,
                fitted = fit)
    class(out) <- "respthresh"
    out
}

### method functions

#' @importFrom stats formula
#' @export
`print.respthresh` <-
    function(x, ...)
{
    cat("Binary threshold for model\n")
    print(formula(x))
    cat("\nthreshold", x$coefficients, "\n")
    cat("average fitted values", x$values, "\n")
    cat("explained deviance", x$expl.deviance, "\n")
    invisible(x)
}

#' @export
`plot.respthresh` <-
    function(x, ...)
{
    plot(x$cutoff, x$devprofile, xlab = "Response Cutoff",
         ylab = "Explained Deviance", type = "l", ...)
    abline(v = x$coefficients, col = 2, ...)
}

#' @importFrom stats deviance qchisq formula
#' @export
`summary.respthresh` <-
    function(object, ...)
{
    tmp <- c(object$null.deviance, deviance(object), object$orig.deviance)
    devtable <- cbind(tmp, c(NA, -diff(tmp)))
    dimnames(devtable) <- list(c("Null", "Threshold", "Model"),
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

#' @importFrom stats printCoefmat
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
    invisible(x)
}

