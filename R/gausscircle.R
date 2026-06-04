### Fit circular Gaussian responses on 2D ordination

#' Circular Gaussian Response on 2D Ordination

#' @param ord Ordination result.
#' @param comm Community data.
#' @param freqlim Frequency limit.
#' @param family Error family.
#' @param unit Force Gaussian responses to unit tolerances.
#' @param choices Ordination axes.
#' @param display Ordination scores.
#' @param \dots Other arguments passed to \code{\link{scores}}.
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
`gaussellipse` <-
    function(ord, comm, freqlim = 10, family = quasipoisson(),
             choices = 1:2, display = "sites", species, ...)
{
    x <- scores(ord, choices = choices, display = display, ...)
    fr <- colSums(comm > 0)
    w <- if (is.atomic(ord)) NULL else weights(ord)
    if (is.null(w)) w <- rep.int(1, nrow(x))
    x <- cbind(x, x^2, x[,1] * x[,2])
    ## use names of Oksanen et al. Ecology 82 (2001), p. 1193, eq. 10
    ## mu = exp(a + b1*x + b2*x^2 + c1*y + c2*y^2 + d*x*y)
    colnames(x) <- c("b1", "c1", "b2", "c2", "d")
    x <- model.matrix( ~ ., as.data.frame(x))
    y <- comm[,species]
    mod <- glm.fit(x, y, family = family, weights = w)
    p <- coef(mod)
    ## new constants
    q0 <- 4 * p["b2"]*p["c2"] - p["d"]^2
    if (q0 <= 0 || p["b2"] + p["c2"] >= 0)
        stop("not a Gaussian surface")
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
    ## return
    out <- c(xopt, yopt, xtol, ytol, rxy)
    names(out) <- c("xopt", "yopt", "xtol", "ytol", "rxy")
    out
}
