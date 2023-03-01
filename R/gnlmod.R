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
    ## error distribution
    Dev <- family()$dev.resids
    V <- family()$variance
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
    out$y <- y
    mu <- eval(SSmodel, envir = split(out$estimate, pnames))
    out$fitted.values <- mu
    out$residuals <- (y - mu) / mu
    out$data <- data
    out$deviance <- 2 * out$minimum
    out
}
