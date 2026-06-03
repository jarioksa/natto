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
### diagonal.
###
### Simple way of adding these is:
###
### tol <- tolerance(ord) # notes on scaling below
### sco <- scores(ord, dis="sp")
### plot(ord)
### for (i in 1:nrow(sco))
###     lines(vegan:::veganCovEllipse(
###         cov = diag(tol[i,1:2]),
###         center = sco[i,1:2]),
###     col=2)

### It seems that in cca/ca we get average 1 species tolerances with
### options scaling = "sites", hill = TRUE (or scaling = -1). In
### decorana() the scaling cannot be changed in scores()/plot(), but
### it should be OK and similar to above. It is essential to use same
### scaling in scores(..., display="species") and tolerance(...,
### which="species").
