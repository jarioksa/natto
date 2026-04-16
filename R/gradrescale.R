### Rescale observed ecological gradient similarly as a decorana axis
#' @title Rescale Ecological Gradient to Constant Community Heterogeneity
#'
#' @details
#'
#' Observed environmental gradient is rescaled in similar way as
#' species axis in \code{\link{rdecorana}} and
#' \code{\link[vegan]{decorana}}. The goal is that weighted averages
#' of species have equal dispersion along the rescaled gradient, but
#' retain the order of sites along environmental gradient, unlike in
#' \code{[r]decorana}.
#'
#' Axis rescaling in \code{[r]decorana} is often poorly
#' understood. Although rescaling is based on site (sample) scores, it
#' is primarily performed on species scores, and site (sample) scores
#' are finally found as weighted averages of rescaled species
#' scores. With this the rank order of site (sample) scores can
#' change. \code{[r]decorana} usually performs several sweeps (default
#' is four), and when new species scores are based on rescaled site
#' (sample) scores, their rank order can also change from the
#' original. Changing the order is not acceptable in rescaling observed
#' gradients, and therefore input gradient is handled like species
#' scores in \code{[r]decorana} and only one rescaling sweep is made.
#'
#' @seealso \code{\link{rdecorana}}, \code{\link[vegan]{decorana}}.
#'
#' @param grad Observed numeric gradient.
#' @param comm Community data matrix.
#'
#' @return Rescaled gradient scaled in units of heterogeneity.
#'
#' @author Jari Oksanen

#' @importFrom vegan wascores
#'
#' @examples
#' if(require(vegan, quietly=TRUE)) {
#' data(mite, mite.env)
#' WatrCont.rescaled <- with(mite.env, gradrescale(WatrCont, mite))
#' plot(WatrCont.rescaled ~ WatrCont, data=mite.env)
#' abline(lm(WatrCont.rescaled ~ WatrCont, data=mite.env), col=4)
#' legend("topleft", c("Nonlinear rescaling", "Linear scaling"),
#'     lty=c(NA,1), pch=c(1,NA), col=c(1,4))
#' }
#'
#' @export
`gradrescale` <-
    function(grad, comm)
{
    wa <- wascores(grad, comm, expand = FALSE)
    comm <- t(comm)
    ## stretch from decorana
    out <- stretch(comm, wa, grad, rowSums(comm))
    grad <- out$cproj
    grad <- grad - min(grad)
    ## grad length was assessed for rproj
    grad <- grad / max(grad) * diff(range(out$rproj))
    grad
}
