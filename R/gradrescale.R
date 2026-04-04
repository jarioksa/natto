### Rescale observed ecological gradient similarly as a decorana axis
#' @title Rescale Ecological Gradient to Constant Community Heterogeneity
#'
#' @details
#' Observed environmental gradient is rescaled similarly as axis in
#' \code{\link{rdecorana}} ahd \code{\link[vegan]{decorana}}. The goal
#' is that weighted averages of species have equal dispersion along
#' the rescaled gradient.
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
