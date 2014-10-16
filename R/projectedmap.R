#' Return a projected map with given extent
#'
#' The function takes either limits of extent or a two-column matrix
#' and returns a projected map containing the original limits or the
#' points.
#'
#' @param extent Either (1) a four element numeric vector of east min,
#' east max, north min, north max or (2) a two-column matrix of
#' coordinates of points that should be shown on the map, or (3) a
#' species occurence data from \code{\link[spocc]{occ}} (\pkg{spocc}
#' package).
#'
#' @param CRS CRS string defining the output projection.
#'
#' @param inCRS CRS string describing the input projection of
#' \code{extent}. This will not be used if CRS is already defined for
#' the \code{extent}.
#'
#' @param resolution Resolution of the map as defined in
#' \code{\link[rworldmap]{getMap}}.
#'
#' @param pad extend climits by this proportion of coordinate range.
#'
#' @import sp
#' @importFrom rgeos gIntersection
#' @importFrom raster extent as.vector
#' @import rgdal
#' @importFrom rworldmap getMap
#' @importFrom spocc occ2df
#'
#' @examples
#' ## Lambert Azimuthal Equal Area: EU recommendation
#' pl <-projectedmap(c(5,31,54,71), CRS="+init=epsg:3035", pad=0)
#' plot(pl, col="gray")
#' ## display extent defined by latitude and longitude ranges
#' llgridlines(pl, east=c(5,31), north=c(54,71), lty=1)
#' 
#' ## Research areas in Virtanen et al., J. Veg. Sci. 17, 519-528 (2006)
#' sites <- structure(list(lon = c(41, 46, 50, 54, 67, 70, 102, 
#' 112, 117, 141, 142, 149, 163, 179, 16, 138), lat =
#' c(67.3, 68.15, 69.15, 68.5, 70.45, 73, 77.2, 76, 73, 72, 75, 71, 69,
#' 71, 78, 75)), .Names = c("longit", "latid"), row.names = c("Kack",
#' "Kani", "Kolg", "Pech", "WYam", "NYam", "Chel",
#' "Neta", "Olon", "Yana", "Fadd", "Lopa", "Koly",
#' "Wran", "Sval", "Kote"), class = "data.frame")
#' ## LAEA North Pole Russia -- default padding will extend the area to an
#' ## illegal polygon
#' pl <- projectedmap(sites, CRS= "+init=epsg:3576", pad=0)
#' plot(pl)
#' points(attr(pl, "input"), pch=16, col=4)

#' #' @export
`projectedmap` <-
    function(extent, CRS = "+proj=longlat +datum=WGS84",
             inCRS = "+proj=longlat +datum=WGS84",
             resolution = "low", pad = 0.04)
{
    NADD <- 21
    ## get map
    map <- getMap(resolution)
    map4 <- proj4string(map)
    ## if extent is a vector of four elements, expand it to two-column
    ## matrix by adding points along the margins so that the map will
    ## completely show defined area.
    if (is.vector(extent) && length(extent) == 4) {
        x <- seq(extent[1], extent[2], length=NADD)
        y <- seq(extent[3], extent[4], length=NADD)
        extent <- data.frame("E" = c(x, rep(extent[2], NADD),
                             rev(x), rep(extent[1], NADD)),
                             "N" = c(rep(extent[3], NADD), y,
                             rep(extent[4], NADD), rev(y)))
    } else if (inherits(extent, "occdat")) {
        ## species occurrence data obtained by spocc::occ() query
        ## belong to class 'occdat'. Change it to a data frame of
        ## longitude and latitude
        extent <- occ2df(extent)[,2:3, drop=FALSE]
    }
    ## now extent should be two-column data.frame or a matrix.
    if (ncol(extent) != 2)
        stop("extent should have two columns or four elements")
    ## take care extent is SpatialPoints object
    if (!inherits(extent, "SpatialPoints")) {
        if (is.matrix(extent))
            extent <- as.data.frame(extent)
            coordinates(extent) <-
                as.formula(paste("~", paste(colnames(extent), collapse = "+")))
    }
    ## use default inCRS if extent has no proj4string
    if (is.na(proj4string(extent)))
        proj4string(extent) <- inCRS
    ## change extent to the output projection, and keep the points
    extent <- spTransform(extent, CRS(CRS))
    ## Look for limits to give the projected extent
    clip <- as.vector(extent(extent))
    ## maps are in equal aspect ratio, and padding is defined in
    ## absolute projected units
    if (pad > 0) {
        pad <- pad * max(diff(clip[1:2]), diff(clip[1:2])) / 2
        clip <- clip + c(-pad, pad, -pad, pad)
    }
    x <- seq(clip[1], clip[2], length=NADD)
    y <- seq(clip[3], clip[4], length=NADD)
    clip <- data.frame("E" = c(x, rep(clip[2], NADD),
                         rev(x), rep(clip[1], NADD)),
                         "N" = c(rep(clip[3], NADD), y,
                         rep(clip[4], NADD), rev(y)))
    coordinates(clip) <- ~ E + N
    proj4string(clip) <- CRS
    ## transform projected rectangle to map coordinates
    clip <- spTransform(clip, CRS(map4))
    ## make clip into SpatialPolygons
    clip <- Polygons(list(Polygon(clip, hole=FALSE)), ID="clip")
    clip <- SpatialPolygons(list(clip), proj4string=CRS(map4))
    map <- gIntersection(map, clip, byid = TRUE)
    if(is.na(proj4string(map)))
       proj4string(map) <- map4
    map <- spTransform(map, CRS(CRS))
    attr(map, "input") <- extent
    map
}
