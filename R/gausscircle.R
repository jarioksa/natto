### Fit Gaussian responses on 2D ordination to estimate species
### tolerances

#' Gaussian Response Circles or Ellipses on 2D Ordination
#'
#' Correspondence Analysis (CA, CCA) is often represented as unimodal
#' method approximating Gaussian responses. Axes of sites can be seen
#' as proxies of environmental gradients, and species scores as
#' estimates of species optima along those gradients. Ter Braak &
#' Looman (1986) demonstrated this model to work well when species
#' have Gaussian responses with equal response widths, known as
#' tolerances. It is possible to scale axes by species tolerance
#' units, so that all species have tolerance 1 on average. Functions
#' \code{gausscircle} and \code{gaussellipse} help in illustrating and
#' inspecting these conjectures by actually fitting Gaussian responses
#' for species and returning their parameters, or the locations of
#' optima and widths of tolerances for each species on 2D ordination
#' graph. Function \code{gausscircle} assumes equal tolerances on both
#' axes (and optionally forces this tolerance to 1), and
#' \code{gaussellipse} fits Gaussian response with independent
#' tolerances and with interaction terms.
#'
#' The functions were written for curiosity, and they do not (yet)
#' have good support functions. If you use these, you must be prepared
#' to learn how to apply the results (see Examples).
#'
#' For the supposed unit-tolerance you must use adequate scaling of
#' ordination. For \code{\link[vegan]{cca}} and
#' \code{\link[vegan]{ca}}, \code{scaling = "sites"} with \code{hill =
#' TRUE} (or numerical shortcut \code{scaling = -1}) should give
#' average tolerance 1, and \code{\link[vegan]{decorana}} scaling
#' should do so automatically. In addition, with rescaling (default)
#' \code{decorana} tries to make the average tolerance 1 all along the
#' axis.
#'
#' \pkg{vegan} function \code{\link[vegan]{tolerance}} for \code{cca}
#' and \code{decorana} uses weighted averages methods to find the
#' tolerances of species, and these can be used and displayed in the
#' same way as the results of these functions.
#'
#' The methods are based on curve fitting and need sufficient non-zero
#' data. The parameters are \code{NA} for species below
#' \code{freqlim}. They are based on fitting Gaussian model as
#' polynomial regression, and translating polynomial coefficients to
#' Gaussian parameters (ter Braak & Looman 1986). The translation in
#' \code{gaussellipse} follows Oksanen et al. (2001). The polynomial
#' models may find a model with no optimum type response, and these
#' results are skipped and returned as \code{NA}.
#'
#' @references
#'
#' Oksanen, J., Läärä, E., Tolonen, K. & Warner, B.G. 2001. Confidence
#' intervals for the optimum in the Gaussian response
#' function. \emph{Ecology} 82, 1191--1197.
#'
#' ter Braak, C.J.F & Looman, C.W.N 1986. Weighted averaging, logistic
#' regression and the Gaussian response model. \emph{Vegetatio} 65,
#' 3--11.

#' @param ord Ordination result. The unimodality is suggested for the
#'     Correspondendence analysis of methods, but the function is
#'     ignorant, and can also analyse other models, even when this
#'     makes no sense (but \code{\link[vegan]{metaMDS}} results with
#'     WA scores for species is a legitimate object).
#' @param comm Community data.
#' @param freqlim Frequency limit, species below this limit are
#'     skipped.
#' @param family Error family.
#' @param unit Fit Gaussian responses of unit tolerances.
#' @param choices Ordination axes.
#' @param display Ordination scores on which responses are fitted. The
#'     Gaussian models are for species, so this should be
#'     \code{"sites"} or in \code{cca} alternatively \code{"lc"}.
#' @param \dots Other arguments passed to \code{\link{scores}}.
#'
#' @return A matrix of Gaussian parameters. \code{xopt} and
#'     \code{yopt} are the estimated location of the species optimum
#'     (ordination scores should approximate these ), \code{tol} the
#'     estimated isometric tolerance (\code{gausscircle}) or
#'     \code{xtol}, \code{ytol} and \code{rxy} tolerances for axes and
#'     their correlation (\code{gaussellipse}), and \code{top} the
#'     estimated height of the response at the optimum (only returned
#'     in \code{gausscircle}), and finally a column for
#'     frequency. Each row is for one species. If Gaussian parameters
#'     are \code{NA}, the species was below \code{freqlim} or the
#'     fitted model was not of Gaussian optimum type.
#'
#' @seealso \code{\link[vegan]{tolerance}} for weighted averages
#'     estimates of tolerance, \code{\link[vegan]{wascores}} for
#'     direct estimation of tolerances with observed gradients, and
#'     \CRANpkg{analogue} package for extensive use of tolerances in
#'     environmental calibration. The response fitting uses
#'     \code{\link{glm}}.
#'
#' @examples
#' library(vegan)
#' data(dune, package = "vegan")
#' ord <- cca(dune)
#' tol0 <- vegan::tolerance(ord, scaling = -1)
#' summary(tol0)
#' tol1 <- gausscircle(ord, dune, scaling = -1, freqlim=6)
#' summary(tol1)
#' tol2 <- gaussellipse(ord, dune, scaling = -1, freqlim=6)
#' summary(tol2)
#' ## plotting is harder, but let us start with the easier case with
#' ## tolerance circles
#' plot(ord, display = "species", scaling = -1)
#' vegan::ordilabel(tol1)  # optima, should be close species score
#' symbols(tol1, circles = tol1[,"tol"], inches = FALSE, add = TRUE)
#' title(main = "Tolerance Circles")
#' ## For other we need to add covariance ellipses. vegan::tolerance
#' ## is simpler because there we have diagonal variance matrix. NB.,
#' ## you must use squared tolerances for covariance ellipses.
#' plot(ord, display = "species", scaling = -1)
#' sco <- scores(ord, display = "species", choices = 1:2, scaling = -1)
#' for (i in 1:nrow(tol0))
#'     lines(vegan:::veganCovEllipse(
#'         cov = diag(tol0[i, 1:2]^2, nrow = 2),
#'         center = sco[i, 1:2]))
#' title(main = "WA Tolerances")
#' ## Ellipses need also off-diagonal component of covariances
#' plot(ord, display = "species", scaling = -1)
#' for (i in 1:nrow(tol2)) {
#'     cv <- diag(tol2[i, 3:4]^2, nrow = 2)
#'     cv[2:3] <- tol2[i, 3] * tol2[i, 4] * tol2[i, 5]
#'     if (!anyNA(cv))
#'         lines(vegan:::veganCovEllipse(cv, tol2[i, 1:2]))
#' }
#' vegan::ordilabel(tol2)
#' title(main = "Tolerance Ellipses")
#'
#'
#' @importFrom stats glm.fit model.matrix coef weights quasipoisson
#' @importFrom vegan scores
#'
#' @author Jari Oksanen
#'
#' @export
`gausscircle` <-
    function(ord, comm, freqlim = 10, family = quasipoisson(), unit = FALSE,
             choices = 1:2, display = "sites", ...)
{
    x <- scores(ord, choices = choices, display = display, ...)
    out <- matrix(NA, ncol(comm), 5)
    colnames(out) <- c("xopt", "yopt", "tol", "top", "freq")
    rownames(out) <- colnames(comm)
    fr <- colSums(comm > 0)
    out[,5] <- fr
    w <- if (is.atomic(ord)) NULL else weights(ord)
    if (is.null(w)) w <- rep.int(1, nrow(x))
    if (unit) {
        off <- -0.5 * rowSums(x^2)
    } else {
        x <- cbind(x, x2 = rowSums(x^2))
        off <- rep.int(0, nrow(x))
    }
    x <- model.matrix( ~ ., as.data.frame(x))
    for(i in 1:ncol(comm)) {
        if (fr[i] < freqlim) next
        y <- comm[,i]
        mod <- glm.fit(x, y, family = family, weights = w, offset = off)
        p <- coef(mod)
        if (unit)
            p <- c(p, -0.5)
        if (p[4] >= 0) next
        tol <- sqrt(-1/2/p[4])
        xopt <- -p[2]/2/p[4]
        yopt <- -p[3]/2/p[4]
        top <- mod$family$linkinv(p[1] + p[2]*xopt + p[3]*yopt +
                                  p[4]*(xopt^2 + yopt^2))
        out[i,1:4] <- c(xopt, yopt, tol, top)
    }
    out
}

###################################################

### Notes for development

### the tol-radius circles centred at (xopt, yopt) can be added to an
### ordination graph with

### plot(ord)
### gm <- gausscircle(ord, comm)
### vegan::ordilabel(gm) # optima
### symbols(gm, circles = gm[,"tol"], inches=FALSE, add = TRUE)

### vegan::tolerance functions use reciprocal averaging ideas of
### estimating SD by axis. These can be added to the plot using
### vegan:::veganCovEllipse(). These are ellipses with principal axes
### parallel to the ordination axes, that is, based on covariance
### matrix with off-diagonal 0 (uncorrelated) and tolerances at
### diagonal. NB., veganCovEllipse draws *covariance* ellipses whereas
### `gausscircle` and `gaussellipse` return non-squared tolerances. My
### understanding is that `vegan::tolerance` functions do the same. So
### we need to square tolerances when we plot them using
### `vegan:::veganCovEllipse`. When we use `symbols(..., circles = )`
### then tolerances should be like they are and not squared.
###
### Simple way of adding these is:
###
### tol <- tolerance(ord) # notes on scaling below
### sco <- scores(ord, dis="sp")
### plot(ord)
### for (i in 1:nrow(sco))
###     lines(vegan:::veganCovEllipse(
###         cov = diag(tol[i,1:2]^2),
###         center = sco[i,1:2]),
###     col=2)

### It seems that in cca/ca we get average 1 species tolerances with
### options scaling = "sites", hill = TRUE (or scaling = -1). In
### decorana() the scaling cannot be changed in scores()/plot(), but
### it should be OK and similar to above. It is essential to use same
### scaling in scores(..., display="species") and tolerance(...,
### which="species").

### `gaussellipse` fits full 2D Gaussian response with interaction
### term. The manipulation of polynomial terms is based on Oksanen et
### al. (2001), Ecology 82, 1191-1197. This also returns non-squared
### tolerances for axes with correlation coefficent `rxy` for the
### response. To plot with vegan:::veganCovEllipse, we need a
### covariance matrix:
###
###    gm <- gaussellipse(ord, comm)
###    cv <- diag(gm[i, c("xtol","ytol")]^2, nrow=2)
###    cv[2:3] <- gm[i,"xtol"] * gm[i,"ytol"] * gm[i,"rxy"]
###    lines(vegan:::veganCovEllipse(cv, gm[i,1:2]))

### Still preliminary and for testing only: handles only one species
#' @importFrom stats glm.fit model.matrix coef weights quasipoisson
#' @importFrom vegan scores
#' @rdname gausscircle
#' @export
`gaussellipse` <-
    function(ord, comm, freqlim = 10, family = quasipoisson(),
             choices = 1:2, display = "sites", ...)
{
    x <- scores(ord, choices = choices, display = display, ...)
    fr <- colSums(comm > 0)
    w <- if (is.atomic(ord)) NULL else weights(ord)
    if (is.null(w)) w <- rep.int(1, nrow(x))
    x <- cbind(x, x^2, x[,1] * x[,2])
    ## use names of Oksanen et al. Ecology 82 (2001), p. 1193, eq. 10
    ## mu = exp(a + b1*x + b2*x^2 + c1*y + c2*y^2 + d*x*y)
    colnames(x) <- c("b1", "c1", "b2", "c2", "d")
    out <- matrix(NA, ncol(comm), 6)
    colnames(out) <- c("xopt", "yopt", "xtol", "ytol", "rxy", "freq")
    rownames(out) <- colnames(comm)
    out[,"freq"] <- fr
    x <- model.matrix( ~ ., as.data.frame(x))
    for (i in 1:ncol(comm)) {
        if (fr[i] < freqlim)
            next
        y <- comm[,i]
        mod <- glm.fit(x, y, family = family, weights = w)
        p <- coef(mod)
        ## new constants
        q0 <- 4 * p["b2"]*p["c2"] - p["d"]^2
        if (q0 <= 0 || p["b2"] + p["c2"] >= 0)
            next
        p1 <- p["d"] - 2*p["b2"]*p["c1"]/p["b1"]
        p2 <- p["d"] - 2*p["c2"]*p["b1"]/p["c1"]
        ## joint optimum on (x,y), eqs. 11 & 12 in Oksanen et al.
        xopt <- -p["b1"]/2/p["b2"] * (1 + p1 * p["d"] / q0)
        yopt <- -p["c1"]/2/p["c2"] * (1 + p2 * p["d"] / q0)
        ## tolerances, eqs 13 & 14 in Oksanen et al.
        xtol <- sqrt(-1/2/p["b2"] * (1 + p["d"]^2 / q0))
        ytol <- sqrt(-1/2/p["c2"] * (1 + p["d"]^2 / q0))
        ## interaction term eq. 15 in Oksanen et al.
        rxy <- p["d"] / sqrt(4 * p["b2"]*p["c2"])
        out[i,1:5] <- c(xopt, yopt, xtol, ytol, rxy)
    }
    out
}

### USER INTERFACE

### Algebra seems to be OK: checked by verification test and
### vegan::ordisurf(..., knots = 2, family = quasipoisson,
### scaling=...) gives contours fitting ellipses. Not too
### user-friendly at the moment, but needs hacker mind, in particular
### in plotting which needs vegan:::veganCovEllipse and transforming
### tolerances to squared tolerances for variance-covariance matrix.

### Which class? Discussion also handles vegan::tolerance with classes
### "tolerance.cca" & "tolerance.decorana", both inheriting from
### "tolerance". The basic class is defined in analogue, but there it
### references specifically Weighted Averages, and there is no real
### inheritance from analogue::tolerance -> vegan::tolerance.cca. It
### is attractive to call these here tolerance.* classes, but the lack
### of real inheritance does not make this useful, for instance
### analogue:::print.tolerance() hard-codes header as "Weighted
### Averages Tolerances", and making a generic print.tolerance() here
### would conflict with analogue (and it is prudent to assume these
### may be used together). Naturally, we could keep the output as a
### simple matrix (like now) and let the potential user to figure out
### what to do with the result.

### Method functions?
###
### If we have class, we need a print for cleaner output. We do not
### need it with matrix. Also summary() works meaningfully with matrix.
###
### The only really needed method is plotting
### circles/ellipses. symbols() works nicely with gausscircles(), but
### vegan:::tolerance.cca/decorana and gaussellipses would plot
### covariance ellipses, and in this file I have used unexported
### vegan:::veganCovEllipse() in a loop over rows of matrix. This or
### something similar could need wrapping to a function (with yet
### unknown name) to add ellipses to a graph or to draw a new
### graph. vegan:::tolerance.cca/decorana add ellipses to centred to
### species scores, but species optima change from ordination in
### ordicircle(), ordiellipse(), and these alternative species scores
### should be added as well. Construction of covariance matrix needs
### some care, and it could be helpful to have a vcov() method
### (generic in stats), which could return either vcov matrix for a
### single row (species) or a list of vcov matrices for all
### species. The vcov matrix is diagonal in
### vegan:::tolerance.cca/decorana and in gausscircle(), but
### gaussellipse has covariance elements filled as well. For single
### species we need there something like:
###
### cv <- diag(x[i, 3:4]^2, nrow=2)     # xtol^2, ytol^2
### cv[2:3] <- x[i,3] * x[i,4] * x[i,5] # xtol * ytol * rxy
