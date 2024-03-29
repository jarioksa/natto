#' Generalized Non-Linear Regression Analysis
#'
#' Non-linear regression (\code{\link{nls}}) \code{\link{selfStart}}
#' model generalized to use exponential \code{\link{family}} error
#' distributions.
#'
#' @details
#'
#' Function combines three \R{} functions for fitting non-linear
#' models with exponential family errors:
#' \enumerate{
#'   \item Non-linear regression is defined by
#'      \code{\link{selfStart}} like in \code{\link{nls}}. In addition to
#'      defining the nonlinear function, these also provide initial values and
#'      derivatives of model parameters.
#'   \item \code{\link{family}} functions are used to define the error
#'      distribution allowing use of other than least squares models. The
#'      \code{family} functions provide the log-likelihood function for fitting,
#'      and allows finding the partial derivatives of model parameters for the
#'      log-likelihood function (see McCullagh & Nelder 1989, p. 41).
#'   \item The log-likelihood function (with partial derivatives) is optimized
#'      with \code{\link{nlm}}.
#' }
#' The result is mostly similar to \code{\link{glm}} object, and can be
#' handled with many \code{glm} method functions, except those that assume
#' linear fit.
#'
#' @references
#' McCullagh, P & Nelder, J.A. (1989) \emph{Generalized Linear Models}, 2nd ed.
#'    Chapman & Hall.
#'
#' @examples
#' ## Compare to vegan species-area models
#' require(vegan) || stop("needs vegan")
#' data(sipoo, sipoo.map)
#' S <- specnumber(sipoo)
#' ## Arrhenius model in vegan with least squares
#' nls(S ~ SSarrhenius(area, k, z), sipoo.map)
#' ## Assume Poisson error
#' (mod <- gnlmod(S ~ SSarrhenius(area, k, z), sipoo.map, family=poisson))
#' ## Arrhenius should be identical to Poisson glm with log-link
#' glm(S ~ log(area), sipoo.map, family=poisson) # (Intercept) = log(k)
#' ## some method functions
#' fitted(mod)
#' predict(mod)
#' predict(mod, newdata = list(area = seq(1, 250, len=31)))
#' summary(mod)
#'
#' @param formula Model formula for a \code{\link{selfStart}}
#' \code{\link{nls}} model.
#' @param data Data for formula.
#' @param family Error distribution of the exponential
#' \code{\link{family}}.
#' @param wts Prior weights
#' @param \dots Other parameters passed to \code{\link{nlm}} that
#' performs the actual regression.

#' @return Function returns the \code{\link{nlm}} result object, but
#'     it is amended with \code{\link{glm}} object related to
#'     \code{\link{family}}.
#'
#' @importFrom stats family gaussian getInitial nlm
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
    out$fitted.values <- as.vector(mu)
    out$residuals <- (y - mu) / mu
    out$prior.weights = wts
    n.ok <- length(y) - sum(wts == 0)
    out$df.null <- n.ok - 1
    out$df.residual <- n.ok - length(out$estimate)
    intercept <- sum(wts * y)/sum(wts)
    out$null.deviance <- sum(fam$dev.resids(y, intercept, wts))
    out$rank <- pnames
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
#' @param object \code{gnlmod} result object.
#' @param newdata New independent data for prediction.
#' @param type Type of predicted values; only \code{"response"} implemented.
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
    function(object, ...)
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

#' @importFrom stats cov2cor
#' @rdname gnlmod
#' @param dispersion The dispersion parameter for the family
#'     used. Either a single numerical value of \code{NULL} (the
#'     default), when it is inferred from the \code{object}.
#' @param correlation logical; if \code{TRUE}, the correlation matrix
#'     of the estimated parameters is returned and printed.
#' @importFrom stats naprint pnorm
#' @export
`summary.gnlmod` <-
    function(object, dispersion = NULL, correlation = FALSE, ...)
{
    df.r <- object$df.residual
    df.f <- object$rank
    if (is.null(dispersion))
        dispersion <-
            if (object$family$family %in% c("poisson", "binomial")) 1
            else if (df.r > 0)
                sum(object$weights * object$residuals^2) / df.r
            else
                NaN

    ## covariance matrix from the hessian
    cf <- object$estimate
    covmat.unscaled <- solve(object$hessian)
    dimnames(covmat.unscaled) <- list(names(cf), names(cf))
    covmat <- dispersion * solve(object$hessian)
    var.cf <- diag(covmat)
    s.err <- sqrt(var.cf)
    zvalue <- cf/s.err
    pvalue <- 2 * pnorm(-abs(zvalue))
    cf.table <- cbind(cf, s.err, zvalue, pvalue)
    dimnames(cf.table) <-
        list(names(cf),
             c("Estimate", "Std. Error", "z value","Pr(>|z|)"))
    ## answer
    keep <- match(c("call","family","deviance", "aic", "df.residual",
                    "null.deviance","df.null", "iterations", "code"),
                  names(object), 0L)
    ans <- c(object[keep],
             list(deviance.resid = residuals(object, type = "deviance"),
		  coefficients = cf.table,
		  dispersion = dispersion,
		  df = c(object$rank, df.r, df.f),
                  cov.unscaled = covmat.unscaled,
		  cov.scaled = covmat))
    if (correlation)
        ans$correlation <- cov2cor(covmat.unscaled)
    class(ans) <- c("summary.gnlmod", "summary.glm")
    ans
}
#' @importFrom stats naprint symnum
#' @export
`print.summary.gnlmod` <-
        function (x, digits = max(3L, getOption("digits") - 3L),
	      symbolic.cor = x$symbolic.cor,
	      signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nCall:\n",
	paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")

    cat("\nCoefficients:\n")
    coefs <- x$coefficients
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                     na.print = "NA", ...)

    ##
    cat("\n(Dispersion parameter for ", x$family$family,
	" family taken to be ", format(x$dispersion), ")\n\n",
	apply(cbind(paste(format(c("Null","Residual"), justify="right"),
                          "deviance:"),
		    format(unlist(x[c("null.deviance","deviance")]),
			   digits = max(5L, digits + 1L)), " on",
		    format(unlist(x[c("df.null","df.residual")])),
		    " degrees of freedom\n"),
	      1L, paste, collapse = " "), sep = "")
    if(nzchar(mess <- naprint(x$na.action))) cat("  (", mess, ")\n", sep = "")
    cat("AIC: ", format(x$aic, digits = max(4L, digits + 1L)),"\n\n",
	"Number of nlm iterations: ", x$iterations, " (nlm code: ", x$code, ")",
	"\n", sep = "")

    correl <- x$correlation
    if(!is.null(correl)) {
# looks most sensible not to give NAs for undefined coefficients
#         if(!is.null(aliased) && any(aliased)) {
#             nc <- length(aliased)
#             correl <- matrix(NA, nc, nc, dimnames = list(cn, cn))
#             correl[!aliased, !aliased] <- x$correl
#         }
	p <- NCOL(correl)
	if(p > 1) {
	    cat("\nCorrelation of Coefficients:\n")
	    if(is.logical(symbolic.cor) && symbolic.cor) {# NULL < 1.7.0 objects
		print(symnum(correl, abbr.colnames = NULL))
	    } else {
		correl <- format(round(correl, 2L), nsmall = 2L,
                                 digits = digits)
		correl[!lower.tri(correl)] <- ""
		print(correl[-1, -p, drop=FALSE], quote = FALSE)
	    }
	}
    }
    cat("\n")
    invisible(x)
}
