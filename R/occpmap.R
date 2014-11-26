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
#' @param extent A four element numeric vector of east min, east max,
#' north min, north max, or other structure accepted by
#' \code{\link{projectedmap}}.
#'
#' @param limit Maximum number of occurrences returned as defined in
#' \code{\link[spocc]{occ}}.
#'
#' @param pruning Remove duplicated points rounded to precision in
#' latitude and longitude (e.g., \code{pruning = 0.1} keeps only cases
#' that are not duplicated with 0.1 degree attitude).
#'
#' @param clip (logical) Clip map to the occurrences of species
#' instead of using the original \code{extent}.
#'
#' @param \dots Arguments passed to \code{\link[spocc]{occ}} and
#' \code{\link{projectedmap}}
#'
#' @importFrom spocc occ occ2df
#' 
#' @export
`occpmap` <- function(query, CRS = "+proj=longlat", extent,
                      limit = 1000, pruning, clip = TRUE, ...)
{
    ## Get map for the occurrence data
    map <- projectedmap(extent, CRS = CRS, ...)
    ## polygon corresponding to the map borders
    geom <- projectedrect2wkt(map)
    ## get occurrence data for the map area
    occs <- occ(query, geometry = geom, limit = limit, ...)
    odf <- occ2df(occs)
    if (!missing(pruning)) {
        dup <- duplicated(round(odf[,2:3]/pruning))
        odf <- odf[!dup,, drop=FALSE]
    }
    ## projected map
    map2 <- projectedmap(odf[,2:3], CRS = CRS, ...)
    if (clip) {
        map <- map2
        attr(map, "input") <- attr(map2, "input")
    }
    plot(map, ...)
    points(attr(map, "input"), ...)
    invisible(map)
}
