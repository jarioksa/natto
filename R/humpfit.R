#' No-interaction Model for Hump-backed Species Richness vs. Biomass
#'
#' Function \code{humpfit} fits a no-interaction model for species
#' richness vs. biomass data (Oksanen 1996). This is a null model that
#' produces a hump-backed response as an artifact of plant size and
#' density.
#'
#' @param mass Biomass.
#' @param spno Species richness.
#' @param family Family of error distribution. Any
#'     \code{\link{family}} can be used, but the link function is
#'     always Fisher's diversity model, and other \code{link}
#'     functions are silently ignored.
#' @param start Vector of starting values for all three parameters.
#'
#' @details
#' The no-interaction model assumes that the humped species richness
#' pattern along biomass gradient is an artifact of plant size and
#' density (Oksanen 1996). For low-biomass sites, it assumes that
#' plants have a fixed size, and biomass increases with increasing
#' number of plants. When the sites becomes crowded, the number of
#' plants and species richness reaches the maximum. Higher biomass is
#' reached by increasing the plant size, and then the number of plants
#' and species richness will decrease. At biomasses below the hump,
#' plant number and biomass are linearly related, and above the hump,
#' plant number is proportional to inverse squared biomass. The number
#' of plants is related to the number of species by the relationship
#' (\code{link} function) from Fisher's log-series (Fisher et
#' al. 1943).
#'
#'  The parameters of the model are:
#'  \enumerate{
#'    \item \code{hump}: the location of the hump on the biomass gradient.
#'    \item \code{scale}: an arbitrary multiplier to translate the biomass
#'    into virtual number of plants.
#'    \item \code{alpha}: Fisher's \eqn{\alpha}{alpha} to translate the
#'    virtual number of plants into number of species.
#'  }
#' The parameters \code{scale} and \code{alpha} are intermingled and
#' this function should not be used for estimating Fisher's
#' \eqn{\alpha}{alpha}.  Probably the only meaningful and interesting
#' parameter is the location of the \code{hump}.
#'
#' Function may be very difficult to fit and easily gets trapped into
#' local solutions, or fails with non-Poisson families, and function
#' \code{profile} should be used to inspect the fitted models. You can
#' use \code{plot.profile}, \code{pairs.profile} for graphical
#' inspection of the profiles, and \code{confint} for the profile
#' based confidence intervals.
#'
#' The original model intended to show that there is no need to
#' speculate about \dQuote{competition} and \dQuote{stress} (Al-Mufti
#' et al. 1977), but humped response can be produced as an artifact of
#' using fixed plot size for varying plant sizes and densities.
#'
#' @return The function returns an object of class \code{"humpfit"}
#' inheriting from class \code{"glm"}. The result object has specific
#' \code{summary}, \code{predict}, \code{plot}, \code{points} and
#' \code{lines} methods. In addition, it can be accessed by the
#' following methods for \code{\link{glm}} objects: \code{\link{AIC}},
#' \code{\link{extractAIC}}, \code{\link{deviance}},
#' \code{\link{coef}}, \code{\link{residuals.glm}} (except \code{type
#' = "partial"}), \code{\link{fitted}}, and perhaps some others.  In
#' addition, function \code{\link[ellipse]{ellipse.glm}} (package
#' \pkg{ellipse}) can be used to draw approximate confidence ellipses
#' for pairs of parameters, if the normal assumptions look
#' appropriate.
#'
#' @author Jari Oksanen

#' @references
#'
#' Al-Mufti, M.M., Sykes, C.L, Furness, S.B., Grime, J.P &
#' Band, S.R. (1977) A quantitative analysis of shoot phenology and
#' dominance in herbaceous vegetation. \emph{Journal of Ecology}
#' 65,759--791.
#'
#' Fisher, R.A., Corbet, A.S. & Williams, C.B. (1943) The relation
#' between the number of species and the number of individuals in a
#' random sample of of an animal population. \emph{Journal of Animal
#' Ecology} 12, 42--58.
#'
#' Oksanen, J. (1996) Is the humped relationship between species
#' richness and biomass an artefact due to plot size? \emph{Journal of
#' Ecology} 84, 293--295.
#'
#' @note
#'
#' The function is a replacement for the original \code{GLIM4}
#' function at the archive of Journal of Ecology.  There the function
#' was represented as a mixed \code{\link{glm}} with one non-linear
#' parameter (\code{hump}) and a special one-parameter link function
#' from Fisher's log-series.  The current function directly applies
#' non-linear maximum likelihood fitting using function
#' \code{\link{nlm}}.  Some expected problems with the current
#' approach are:
#'
#' \itemize{
#'    \item The function is discontinuous at \code{hump} and may be
#'    difficult to optimize in some cases (the lines will always join, but
#'    the derivative jumps).
#'    \item The function does not try very hard to find sensible starting
#'    values and can fail. The user may supply starting values in
#'    argument \code{start} if fitting fails.
#'    \item The estimation is unconstrained, but both \code{scale} and
#'    \code{alpha} should always be positive.  Perhaps they should be
#'    fitted as logarithmic. Fitting \code{\link{Gamma}}
#'    \code{\link{family}} models might become easier, too.
#'    }
#'
#'
#' @seealso \code{\link[vegan]{fisherfit}}, \code{\link[stats]{profile.glm}},
#'     \code{\link[stats]{confint.glm}}.
#'
#' @examples
#'
#' ## Data of Al-Mufti et al. (1977)
#' mass <- c(140,230,310,310,400,510,610,670,860,900,1050,1160,1900,2480)
#' spno <- c(1,  4,  3,  9, 18, 30, 20, 14,  3,  2,  3,  2,  5,  2)
#' (sol <- humpfit(mass, spno))
#' summary(sol) # Almost infinite alpha...
#' plot(sol)
#' ## Parameter 3 (alpha) is too extreme for profile and confint, and we
#' ## must use only "hump" and "scale".
#' plot(profile(sol, parm=1:2))
#' confint(sol, parm=c(1,2))

#'
#' @importFrom stats family nlm poisson
#' @rdname humpfit

#' @export
`humpfit` <-
    function (mass, spno, family = poisson, start)
{
    hump <- function(p, mass, spno, ...) {
        x <- ifelse(mass < p[1], mass/p[1], p[1] * p[1]/mass/mass)
        fv <- p[3] * log(1 + p[2] * x/p[3])
        n <- wt <- rep(1, length(x))
        dev <- sum(dev.resids(spno, fv, wt))
        aicfun(spno, n, fv, wt, dev)/2
    }
    fam <- family(link = "identity")
    aicfun <- fam$aic
    dev.resids <- fam$dev.resids
    if (missing(start))
        p <- c(mean(mass), 100, 10)
    else
        p <- start
    fit <- nlm(hump, p = p, mass = mass, spno = spno, hessian = TRUE)
    p <- fit$estimate
    names(p) <- c("hump", "scale", "alpha")
    x <- ifelse(mass < p[1], mass/p[1], p[1] * p[1]/mass/mass)
    fv <- p[3] * log(1 + p[2] * x/p[3])
    res <- dev.resids(spno, fv, rep(1, length(x)))
    dev <- sum(res)
    residuals <- spno - fv
    aic <- fit$minimum * 2 + 6
    rdf <- length(x) - 3
    out <- list(nlm = fit, family = fam, y = spno, x = mass,
                coefficients = p, fitted.values = fv, aic = aic, rank = 3,
                df.residual = rdf, deviance = dev, residuals = residuals,
                prior.weights = rep(1, length(x)))
    class(out) <- c("humpfit", "glm")
    out
}
#' @param segments Number of segments used for lines.
#' @importFrom stats fitted predict
#' @importFrom graphics lines
#' @rdname humpfit
#' @export
`lines.humpfit` <-
    function(x, segments=101,  ...)
{
    mass <- x$x
    if (!is.null(segments) && segments > 0) {
        mass <- seq(min(mass), max(mass), length=segments)
        fv <- predict(x, newdata = mass)
    }
    else {
        i <- order(mass)
        fv <- fitted(x)
        mass <- mass[i]
        fv <- fv[i]
    }
    lines(mass, fv, ...)
    invisible()
}
#' @param xlab,ylab Axis labels.
#' @param lwd Line width.
#' @param l.col,p.col Line and point colours.
#' @param type Type of ype of \code{plot}: \code{"p"} for observed
#'     points, \code{"l"} for fitted lines, \code{"b"} for both, and
#'     \code{"n"} for only setting axes.
#' @importFrom graphics plot points lines
#' @rdname humpfit
#' @export
`plot.humpfit` <-
    function(x, xlab="Biomass", ylab="Species Richness", lwd=2, l.col="blue",
             p.col = 1, type="b", ...)
{
    plot(x$x, x$y, xlab = xlab, ylab = ylab, type="n", ...)
    if (type == "b" || type == "p")
        points(x, col = p.col, ...)
    if (type == "b" || type == "l")
        lines(x, lwd = lwd, col = l.col, ...)
    invisible()
}
#' @param x Fitted result object.
#' @param \dots Other parameters to functions.
#' @importFrom graphics points
#' @rdname humpfit
#' @export
`points.humpfit` <-
    function(x, ...)
{
    points(x$x, x$y, ...)
    invisible()
}
#' @param newdata Values of \code{mass} used in \code{predict}. The
#'     original data values are used if missing.
#' @importFrom stats fitted coef
#' @rdname humpfit
#' @export
`predict.humpfit` <-
    function(object, newdata = NULL, ...)
{
    if (is.null(newdata))
        return(fitted(object))
    else {
        p <- coef(object)
        x <- unlist(newdata)
        x <- ifelse(x < p[1], x/p[1], p[1]*p[1]/x/x)
        fv <- p[3]*log(1 + p[2]*x/p[3])
    }
    fv
}
#' @importFrom stats family coef deviance df.residual
#' @export
`print.humpfit` <-
    function(x, ...)
{
    cat("\nHump-backed Null model of richness vs. productivity\n\n")
    cat("Family:", family(x)$family,"\n")
    cat("Link function: Fisher diversity\n\n")
    cat("Coefficients:\n\n")
    print(coef(x))
    cat("\nDeviance", deviance(x), "with", df.residual(x))
    cat(" residual degrees of freedom\n")
    invisible(x)
}
#' @importFrom stats printCoefmat
#' @export
`print.summary.humpfit` <-
    function (x, ...)
{
    cat("\nHump-backed Null model of richness vs. productivity\n\n")
    cat("Family:", x$family, "\n")
    cat("Link function: Fisher diversity\n\n")
    cat("Coefficients:\n\n")
    printCoefmat(x$est, ...)
    cat("\nDispersion parameter for", x$family, "family taken to be", x$dispersion,"\n")
    cat("\nDeviance", x$deviance, "with", x$df.residual)
    cat(" residual degrees of freedom\n")
    cat("AIC:", x$aic, "   BIC:", x$bic, "\n")
    cat("\nCorrelation of Coefficients:\n")
    correl <- format(round(x$correlation, 2), nsmall = 2)
    correl[!lower.tri(correl)] <- ""
    print(correl[-1, -3], quote = FALSE)
    cat("\nDiagnostics from nlm:\n")
    cat("Number of iterations: ", x$iter, ", code: ", x$code,
        "\n", sep = "")
    invisible(x)
}
#' @param fitted Fitted result object.
#' @param parm Profiled parameters.
#' @param alpha,maxsteps,del Parameters for profiling range and
#'     density (see Details).
#' @importFrom stats coefficients qchisq qf nlm
#' @rdname humpfit
#' @export
`profile.humpfit` <-
    function(fitted, parm=1:3, alpha=0.01, maxsteps = 20, del = zmax/5, ...)
{
    INSERT3 <- function(vec, fix, val) {
        switch(fix,
               c(val, vec),
               c(vec[1], val, vec[2]),
               c(vec, val)
               )
    }
    HUMP <- function(p, mass, spno, fix, val, ...) {
        b <- INSERT3(p, fix, val)
        x <- ifelse(mass < b[1], mass/b[1], b[1]*b[1]/mass/mass)
        fv <- b[3] * log(1 + b[2]*x/b[3])
        n <- wt <- rep(1, length(x))
        dev <- sum(dev.resids(spno, fv, wt))
        aicfun(spno, n, fv, wt, dev)/2
    }
    dev.resids <- fitted$family$dev.resids
    aicfun <- fitted$family$aic
    minll <- fitted$nlm$minimum
    p <- coefficients(fitted)
    pv0 <- t(as.matrix(p))
    n <- length(fitted$residuals)
    Pnames <- names(p)
    summ <- summary(fitted)
    dispersion <- summ$dispersion
    std.err <- summ$est[, "Std. Error"]
    if (summ$family == "poisson") {
        zmax <- sqrt(qchisq(1 - alpha/2, 3))
        profName <- "z"
    } else {
        zmax <- sqrt(3 * qf(1 - alpha/2, 3, n-3))
        profName <- "tau"
    }
    prof <- vector("list", length = length(parm))
    names(prof) <- Pnames[parm]
    for (i in parm) {
        zi <- 0
        par <- pv0
        pvi <- pv0[-i]
        pi <- Pnames[i]
        for (sgn in c(-1, 1)) {
            step <- 0
            z <- 0
            while ((step <- step + 1) < maxsteps && abs(z) < zmax) {
                bi <- p[i] + sgn * step * del * std.err[i]
                fm <- nlm(HUMP, p = pvi, mass = fitted$x, spno = fitted$y, fix = i, val = bi)
                pvi <- fm$estimate
                ri <- INSERT3(pvi, i, bi)
                names(ri) <- Pnames
                par <- rbind(par, ri)
                zz <- 2*(fm$minimum - minll)/dispersion
                if (zz > -0.001)
                    zz <- max(0, zz)
                else
                    stop(gettextf(
                        "profiling found a better solution:  original fit had not converged:\n%s: %f",
                        Pnames[i], bi))
                z <- sgn*sqrt(zz)
                zi <- c(zi, z)
            }
        }
        si <- order(zi)
        prof[[pi]] <- structure(data.frame(zi[si]), names= profName)
        prof[[pi]]$par.vals <- par[si,]
    }
    val <- structure(prof, original.fit = fitted, summary = summ)
    class(val) <- c("profile.humpfit", "profile.glm", "profile")
    val
}
#' @param object Fitted result object.
#' @importFrom stats coef AIC family deviance df.residual
#' @rdname humpfit
#' @export
`summary.humpfit` <-
    function (object, ...)
{
    dispersion <- if (any(object$family$family == c("binomial",
                          "poisson")))
        1
    else sum(object$residuals^2)/object$df.residual
    p <- coef(object)
    se <- sqrt(dispersion * diag(solve(object$nlm$hessian)))
    est <- cbind(p, se)
    colnames(est) <- c("Estimate", "Std. Error")
    covmat <- solve(object$nlm$hessian)
    dg <- sqrt(diag(covmat))
    cormat <- covmat/outer(dg, dg)
    colnames(cormat) <- names(p)
    rownames(cormat) <- names(p)
    aic <- AIC(object)
    bic <- AIC(object, k = log(length(object$y)))
    out <- list(est = est, aic = aic, bic = bic, family = family(object)$family,
                deviance = deviance(object), df.residual = df.residual(object),
                dispersion = dispersion, correlation = cormat,
                cov.unscaled = covmat,  iter = object$nlm$iterations,
                code = object$nlm$code)
    class(out) <- "summary.humpfit"
    out
}
#' @param segments Number of segments used for lines.
#' @importFrom stats fitted predict
#' @importFrom graphics lines
#' @rdname humpfit
#' @export
`lines.humpfit` <-
    function(x, segments=101,  ...)
{
    mass <- x$x
    if (!is.null(segments) && segments > 0) {
        mass <- seq(min(mass), max(mass), length=segments)
        fv <- predict(x, newdata = mass)
    }
    else {
        i <- order(mass)
        fv <- fitted(x)
        mass <- mass[i]
        fv <- fv[i]
    }
    lines(mass, fv, ...)
    invisible()
}
#' @param xlab,ylab Axis labels.
#' @param lwd Line width.
#' @param l.col,p.col Line and point colours.
#' @param type Type of ype of \code{plot}: \code{"p"} for observed
#'     points, \code{"l"} for fitted lines, \code{"b"} for both, and
#'     \code{"n"} for only setting axes.
#' @importFrom graphics plot points lines
#' @rdname humpfit
#' @export
`plot.humpfit` <-
    function(x, xlab="Biomass", ylab="Species Richness", lwd=2, l.col="blue",
             p.col = 1, type="b", ...)
{
    plot(x$x, x$y, xlab = xlab, ylab = ylab, type="n", ...)
    if (type == "b" || type == "p")
        points(x, col = p.col, ...)
    if (type == "b" || type == "l")
        lines(x, lwd = lwd, col = l.col, ...)
    invisible()
}
#' @param x Fitted result object.
#' @param \dots Other parameters to functions.
#' @importFrom graphics points
#' @rdname humpfit
#' @export
`points.humpfit` <-
    function(x, ...)
{
    points(x$x, x$y, ...)
    invisible()
}
#' @param newdata Values of \code{mass} used in \code{predict}. The
#'     original data values are used if missing.
#' @importFrom stats fitted coef
#' @rdname humpfit
#' @export
`predict.humpfit` <-
    function(object, newdata = NULL, ...)
{
    if (is.null(newdata))
        return(fitted(object))
    else {
        p <- coef(object)
        x <- unlist(newdata)
        x <- ifelse(x < p[1], x/p[1], p[1]*p[1]/x/x)
        fv <- p[3]*log(1 + p[2]*x/p[3])
    }
    fv
}
#' @importFrom stats family coef deviance df.residual
#' @export
`print.humpfit` <-
    function(x, ...)
{
    cat("\nHump-backed Null model of richness vs. productivity\n\n")
    cat("Family:", family(x)$family,"\n")
    cat("Link function: Fisher diversity\n\n")
    cat("Coefficients:\n\n")
    print(coef(x))
    cat("\nDeviance", deviance(x), "with", df.residual(x))
    cat(" residual degrees of freedom\n")
    invisible(x)
}
#' @importFrom stats printCoefmat
#' @export
`print.summary.humpfit` <-
    function (x, ...)
{
    cat("\nHump-backed Null model of richness vs. productivity\n\n")
    cat("Family:", x$family, "\n")
    cat("Link function: Fisher diversity\n\n")
    cat("Coefficients:\n\n")
    printCoefmat(x$est, ...)
    cat("\nDispersion parameter for", x$family, "family taken to be", x$dispersion,"\n")
    cat("\nDeviance", x$deviance, "with", x$df.residual)
    cat(" residual degrees of freedom\n")
    cat("AIC:", x$aic, "   BIC:", x$bic, "\n")
    cat("\nCorrelation of Coefficients:\n")
    correl <- format(round(x$correlation, 2), nsmall = 2)
    correl[!lower.tri(correl)] <- ""
    print(correl[-1, -3], quote = FALSE)
    cat("\nDiagnostics from nlm:\n")
    cat("Number of iterations: ", x$iter, ", code: ", x$code,
        "\n", sep = "")
    invisible(x)
}
#' @param fitted Fitted result object.
#' @param parm Profiled parameters.
#' @param alpha,maxsteps,del Parameters for profiling range and
#'     density (see Details).
#' @importFrom stats coefficients qchisq qf nlm
#' @rdname humpfit
#' @export
`profile.humpfit` <-
    function(fitted, parm=1:3, alpha=0.01, maxsteps = 20, del = zmax/5, ...)
{
    INSERT3 <- function(vec, fix, val) {
        switch(fix,
               c(val, vec),
               c(vec[1], val, vec[2]),
               c(vec, val)
               )
    }
    HUMP <- function(p, mass, spno, fix, val, ...) {
        b <- INSERT3(p, fix, val)
        x <- ifelse(mass < b[1], mass/b[1], b[1]*b[1]/mass/mass)
        fv <- b[3] * log(1 + b[2]*x/b[3])
        n <- wt <- rep(1, length(x))
        dev <- sum(dev.resids(spno, fv, wt))
        aicfun(spno, n, fv, wt, dev)/2
    }
    dev.resids <- fitted$family$dev.resids
    aicfun <- fitted$family$aic
    minll <- fitted$nlm$minimum
    p <- coefficients(fitted)
    pv0 <- t(as.matrix(p))
    n <- length(fitted$residuals)
    Pnames <- names(p)
    summ <- summary(fitted)
    dispersion <- summ$dispersion
    std.err <- summ$est[, "Std. Error"]
    if (summ$family == "poisson") {
        zmax <- sqrt(qchisq(1 - alpha/2, 3))
        profName <- "z"
    } else {
        zmax <- sqrt(3 * qf(1 - alpha/2, 3, n-3))
        profName <- "tau"
    }
    prof <- vector("list", length = length(parm))
    names(prof) <- Pnames[parm]
    for (i in parm) {
        zi <- 0
        par <- pv0
        pvi <- pv0[-i]
        pi <- Pnames[i]
        for (sgn in c(-1, 1)) {
            step <- 0
            z <- 0
            while ((step <- step + 1) < maxsteps && abs(z) < zmax) {
                bi <- p[i] + sgn * step * del * std.err[i]
                fm <- nlm(HUMP, p = pvi, mass = fitted$x, spno = fitted$y, fix = i, val = bi)
                pvi <- fm$estimate
                ri <- INSERT3(pvi, i, bi)
                names(ri) <- Pnames
                par <- rbind(par, ri)
                zz <- 2*(fm$minimum - minll)/dispersion
                if (zz > -0.001)
                    zz <- max(0, zz)
                else
                    stop(gettextf(
                        "profiling found a better solution:  original fit had not converged:\n%s: %f",
                        Pnames[i], bi))
                z <- sgn*sqrt(zz)
                zi <- c(zi, z)
            }
        }
        si <- order(zi)
        prof[[pi]] <- structure(data.frame(zi[si]), names= profName)
        prof[[pi]]$par.vals <- par[si,]
    }
    val <- structure(prof, original.fit = fitted, summary = summ)
    class(val) <- c("profile.humpfit", "profile.glm", "profile")
    val
}
#' @param object Fitted result object.
#' @importFrom stats coef AIC family deviance df.residual
#' @rdname humpfit
#' @export
`summary.humpfit` <-
    function (object, ...)
{
    dispersion <- if (any(object$family$family == c("binomial",
                          "poisson")))
        1
    else sum(object$residuals^2)/object$df.residual
    p <- coef(object)
    se <- sqrt(dispersion * diag(solve(object$nlm$hessian)))
    est <- cbind(p, se)
    colnames(est) <- c("Estimate", "Std. Error")
    covmat <- solve(object$nlm$hessian)
    dg <- sqrt(diag(covmat))
    cormat <- covmat/outer(dg, dg)
    colnames(cormat) <- names(p)
    rownames(cormat) <- names(p)
    aic <- AIC(object)
    bic <- AIC(object, k = log(length(object$y)))
    out <- list(est = est, aic = aic, bic = bic, family = family(object)$family,
                deviance = deviance(object), df.residual = df.residual(object),
                dispersion = dispersion, correlation = cormat,
                cov.unscaled = covmat,  iter = object$nlm$iterations,
                code = object$nlm$code)
    class(out) <- "summary.humpfit"
    out
}
