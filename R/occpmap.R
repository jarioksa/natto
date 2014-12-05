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
#' @param resolution Resolution of the map as defined in
#' \code{\link[rworldmap]{getMap}}.
#'
#' @param border Border colour in the map.
#'
#' @param land Land colour in the map.
#'
#' @param pch Plotting character for points.
#'
#' @param pcol Point colour.
#'
#' @param annotate Add species name as the title and projection as the
#' right margin text in the map.
#' 
#' @param \dots Arguments passed to \code{\link[spocc]{occ}} and
#' \code{\link{projectedmap}}.
#'
#' @importFrom spocc occ occ2df
#' 
#' @export
`occpmap` <- function(query, CRS = "+proj=longlat", extent,
                      limit = 1000, pruning, clip = TRUE, resolution = "low",
                      border = 1, land="grey90", pch=16, pcol = 2,
                      annotate = TRUE,
                      ...)
{
    ## Get map for the occurrence data
    map <- projectedmap(extent, CRS = CRS, resolution = resolution, ...)
    ## polygon corresponding to the map borders
    geom <- projectedrect2wkt(map)
    ## get occurrence data for the map area
    occs <- occ(query, geometry = geom, limit = limit, ...)
    odf <- occ2df(occs)
    taxonname <- unique(odf[,1])
    if (length(taxonname) > 1) {
        warning("records contain several names, using first: ", taxonname)
        taxonname <- taxonname[1]
    }
    if (NROW(odf) == 0)
        stop("no data for ", query)
    if (NROW(odf) == limit)
        warning("number of records equals limit (", limit, "): data may be truncated")
    if (!missing(pruning)) {
        dup <- duplicated(round(odf[,2:3]/pruning))
        odf <- odf[!dup,, drop=FALSE]
    }
    ## projected map
    map2 <- projectedmap(odf[,2:3], CRS = CRS, resolution = resolution, ...)
    if (clip) {
        map <- map2
        attr(map, "input") <- attr(map2, "input")
    }
    plot(map, col = land, border = border, ...)
    points(attr(map, "input"), pch = pch, col = pcol, ...)
    if (annotate) {
        mtext(proj4string(map), 4, line=-3.5, cex=0.8)
        title(main=taxonname, line=-1)
    }
    invisible(map)
}
