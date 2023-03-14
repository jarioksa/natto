#' Generalized Non-Linear Regression Analysis
#'
#' Non-linear regression (\code{\link{nls}}) \code{\link{selfStart}}
#' model generalized to use exponential \code{\link{family}} error
#' distributions
#'
#' @param formula Model formula for a \code{\link{selfStart}}
#' \code{\link{nls}} model.
#' @param data Data for formula.
#' @param family Error distribution of the exponential
#' \code{\link{family}}.
#' @param wts Prior weights
#' @param \dots Other parameters passed to \code{\link{nlm}} that
#' performs the actual regression.

#' @return The \code{\link{nlm}} result object augmented with some
#'     elements of \code{\link{glm}} object.
#'
#' @importFrom stats family getInitial nlm
#'
#' @export
`gnlmod` <-
    function(formula, data, family = gaussian, wts, ...)
{
    if (missing(data))
        data <- NULL
    ## extract dependent data and selfStart model
    y <- eval(formula[[2]], data)
    x <- formula[[3]][[2]]
    assign(as.character(x), eval(x, data))
    SSmodel <- formula[[3]]
    ## starting values
    p <- getInitial(formula, data)
    pnames <- names(p)
    ## family
    fam <- family(link="identity")
    Dev <- fam$dev.resids
    V <- fam$variance
    if (missing(wts))
        wts <- rep(1, length(y))
    ## loss function
    loss <- function(p, y = y, SSmodel = SSmodel, Dev, V, wts = wts,
                     nm = pnames)
    {
        psplit <- split(p, nm)
        mu <- eval(SSmodel, envir=split(p, nm))
        ll <- sum(Dev(y, mu, wts))/2
        attr(ll, "gradient") <-
            -t(wts * (y - mu) / V(mu)) %*% attr(mu, "gradient")
        ll
    }
    out <- nlm(loss, p = p, y = y, SSmodel = SSmodel, Dev, V, wts = wts,
               hessian = TRUE, ...)
    names(out$estimate) <- pnames
    out$y <- y
    out$nobs <- length(y)
    mu <- eval(SSmodel, envir = split(out$estimate, pnames))
    out$fitted.values <- mu
    out$residuals <- (y - mu) / mu
    out$prior.weights = wts
    n.ok <- length(y) - sum(wts == 0)
    out$df.null <- n.ok - 1
    out$df.residual <- n.ok - length(out$estimate)
    intercept <- sum(wts * y)/sum(wts)
    out$null.deviance <- sum(fam$dev.resids(y, intercept, wts))
    out$aic <- fam$aic(y, length(y), mu, wts, 2 * out$minimum) +
        2 * length(p)
    out$family <- fam
    out$data <- data
    out$formula <- formula
    out$deviance <- 2 * out$minimum
    out$call <- match.call()
    class(out) <- c("gnlmod", "glm")
    out
}

##
#' @rdname gnlmod
#' @export
`predict.gnlmod` <-
    function(object, newdata, type = "response", ...)
{
    if (type != "response")
        stop('only type = "response" implemented')
    if (missing(newdata)) return(as.vector(fitted(object)))
    as.vector(eval(object$formula[[3]],
                   as.list(c(newdata,
                             split(object$estimate, names(object$estimate))))))
}

#' @importFrom stats coef
#' @export
`coef.gnlmod` <-
    function(object)
        object$estimate

## print is edited from stats:::print.glm
#' @export
`print.gnlmod` <-
    function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
    cat("\nCall:  ",
	paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    if(length(coef(x))) {
        cat("Coefficients")
        cat(":\n")
        print.default(format(coef(x), digits = digits),
                      print.gap = 2, quote = FALSE)
    } else cat("No coefficients\n\n")
    cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ",
        x$df.residual, "Residual\n")
    if(nzchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep = "")
    cat("Null Deviance:	   ",	format(signif(x$null.deviance, digits)),
	"\nResidual Deviance:", format(signif(x$deviance, digits)),
	"\tAIC:", format(signif(x$aic, digits)))
    cat("\n")
    invisible(x)
}
