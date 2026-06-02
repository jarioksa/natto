### Fit circular Gaussian responses on 2D ordination

#' @param ord Ordination result.
#' @param comm Community data.
#' @param freqlim Frequency limit.
#' @param family Error family.
#' @param unit Force Gaussian responses to unit tolerances.
#' @param choices Ordination axes.
#' @param display Ordination scores.
#' @param \dots Other arguments passed to \code{\link{scores}}.
#'
#' @importFrom stats glm.fit model.matrix coef weights
#' @importFrom vegan scores
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
    w <- weights(ord) %||% rep.int(1, nrow(x))
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
