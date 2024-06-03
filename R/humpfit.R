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
#' (\code{link} afunction) from Fisher's log-series (Fisher et
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
#' based confidence intervals. (With \R{} prior to 4.4-0 you must use
#' \code{library(MASS)} to access these functions.)
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
#' @seealso \code{\link{fisherfit}}, \code{\link[MASS]{profile.glm}},
#'     \code{\link[MASS]{confint.glm}}.
#'
#' @examples
#'
#' mass <- c(140,230,310,310,400,510,610,670,860,900,1050,1160,1900,2480)
#' spno <- c(1,  4,  3,  9, 18, 30, 20, 14,  3,  2,  3,  2,  5,  2)
#' sol <- humpfit(mass, spno)
#' summary(sol) # Almost infinite alpha...
#' plot(sol)
#' ## confint is in MASS, and impicitly calls profile.humpfit.
#' ## Parameter 3 (alpha) is too extreme for profile and confint, and we
#' ## must use only "hump" and "scale".
#' library(MASS)
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
