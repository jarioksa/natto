#' Returns a WKT polygon in long/lat of a projected rectangle
#'
#' Function returns a WKT polygon that gives the limits of a projected
#' map in longitude/latitude. This polygon can be used to force query
#' of species occurrence query to the area of projected map instead of
#' longitude/latitude rectangle (that is not a rectangle in a
#' projected map).
#' 
#' @param object An object with an extent in projected coordinates
#'
#' @param knots Number of knots filled between corner points of the
#' bounding points when defining the polygon.
#'
#' @importFrom raster extent as.vector
#' @import rgdal
#' @import sp
#' 
#' @export
`projectedrect2wkt` <-
    function(object, knots = 7, ...)
{
    ext <- as.vector(extent(object))
    ## A projected rectangle in map unitis
    x <- seq(ext[1], ext[2], length=knots)
    y <- seq(ext[3], ext[4], length=knots)
    xy <- data.frame("E" = c(x, rep(ext[2], knots-1),
                     rev(x), rep(ext[1], knots-1)),
                     "N" = c(rep(ext[3], knots-1), y,
                     rep(ext[4], knots-1), rev(y)))
    coordinates(xy) <- ~ E + N
    proj4string(xy) <- proj4string(object)
    ## project to long/lat
    xy <- spTransform(xy, CRS("+proj=longlat"))
    xy <- as.data.frame(xy)
    ## pack results to a WKT polygon
    xy <-  paste(apply(apply(xy, 1, paste), 2, paste, collapse=" "),
                 collapse=", ")
    paste("POLYGON((", xy, "))")
}
