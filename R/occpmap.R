#' Plot a projected map of species occurrences
#'
#' Function queries species occurrence data with
#' \code{\link[spocc]{occ}} function and draws a projected map with
#' these occurrences using \code{\link{projectedmap}}.
#'
#' @param query Species to be plotted.
#'
#' @param CRS Projection used in \code{\link{projectedmap}}.
#'
#' @param geometry As defined in \code{\link[spocc]{occ}}, typically
#' corner points of the longitude-latitude polygon.
#'
#' @param limit Maximum number of occurrences returned as defined in
#' \code{\link[spocc]{occ}}.
#'
#' @param pruning Remove duplicated points rounded to precision in
#' latitude and longitude (e.g., \code{pruning = 0.1} keeps only cases
#' that are not duplicated with 0.1 degree attitude).
#'
#' @param \dots Arguments passed to \code{\link[spocc]{occ}} and
#' \code{\link{projectedmap}}
#'
#' @importFrom spocc occ occ2df
#' 
#' @export
`occpmap` <- function(query, CRS = "+proj=longlat", geometry = NULL,
                      limit = 1000, pruning, ...)
{
    ## get occurrence data
    occs <- occ(query, geometry = geometry, limit = limit, ...)
    odf <- occ2df(occs)
    if (!missing(pruning)) {
        dup <- round(odf[,2:3]/pruning)
        odf <- odf[!dup,, drop=FALSE]
    }
    ## projected map
    map <- projectedmap(odf[,2:3], CRS = CRS, ...)
    plot(map, ...)
    points(attr(map, "input"), ...)
    invisible(map)
}
