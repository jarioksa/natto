### Polar Ordination A.K.A. Bray-Curtis Ordination according to Bruce
### McCune & Jim Grace (2002), Analysis of Ecological Communities,
### Chapter 17

#' Polar (Bray-Curtis) Ordination
#'
#' Polar or Bray-Curtis ordination is a historic ordination method
#' that could be performed without computers with simple hand
#' calculations (Bray & Curtis 1957). Ordination axis is found by
#' selecting two extreme points and projecting all points between
#' these end points. The current function follows Beals (1984) in
#' selecting the endpoints, projection of points on axis, and defining
#' the residual distances for later axes.

#' @param d Dissimilarities or distances: a \code{\link{dist}} object.
#' @param k Number of dimensions.
#'
#' @details
#'
#' The implementation follows McCune & Grace (2002, Chapter 17). The
#' endpoints are found using the variance-regression method of Beals
#' (1984). The first endpoint has the highest variance of distances to
#' other points. This guarantees that the point is at the margin of
#' the multivariate cloud but is not an outlier, since outliers have
#' long distances to all points and hence low variance of
#' distances. The second endpoint has the lowest (most negative)
#' regression coefficient between distances from the first and second
#' point to all other points. This selects a point at the margin of
#' the main cloud of points, opposite to first endpoint. All points
#' are projected on the axis between the endpoints using Euclidean
#' geometry, and this gives the scores on a polar ordination
#' axis. Then the effect of the  axis are removed by calculating
#' residual distances using Euclidean geometry. Ecological indices are
#' usually semimetric, and negative residual distances can emerge, but
#' these are taken as zero in the current function.
#'
#' Polar ordination is a historical method that is little used today,
#' but several authors claim that it is a powerful method (see McCune
#' & Grace 2002). Although the basic operations can be easily
#' performed by hand or graphically, the later developments of
#' endpoint selection require more extensive calculations. With modern
#' numerical utilities, the polar ordination is not faster than metric
#' multidimensional scaling (\code{\link{cmdscale}},
#' \code{\link[vegan]{wcmdscale}}).
#'
#' It is possible to use predefined endpoints instead of automatic
#' selection. This can be useful for confirmatory analysis (McCune &
#' Grace 2002). However, this is not (yet) implemented in this
#' function (but contributions are welcome).
#'
#' @return
#'
#' The function returns an object of class \code{"polarord"} with the
#' following elements:
#' \itemize{
#'   \item \code{points}: The ordination scores.
#'   \item \code{inertia}: Total inertia of the input dissimilarities.
#'   \item \code{eig}: Eigenvalues of axes. These do not usually add up to
#'      total inertia and may not be in strictly descending order.
#'   \item \code{endpoints}: The indices (not the names) of the endpoints for
#'      each axis.
#' }
#'
#' @references
#'
#' Beals, E. W. (1984) Bray-Curtis ordination: an effective strategy
#' for analysis of ecological multivariate data. \emph{Advances in
#' Ecological Research} 14, 1--55.
#'
#' Bray, J. R. & Curtis, J. T. (1957) An ordination of the upland
#' forest communities in southern Wisconsin. \emph{Ecological
#' Monographs} 27, 325--349.
#'
#' McCune, B. & Grace, J. B. (2002) \emph{Analysis of Ecological
#' Communities.} MjM Software Design.
#'
#' @examples
#'
#' data(spurn)
#' dis <- dist(spurn, method = "binary") ## Jaccard index
#' ord <- polarord(dis)
#' ord
#' summary(eigenvals(ord))
#' ## add species scores
#' sppscores(ord) <- spurn
#' plot(ord)
#'
#' @importFrom stats var cov dist
#'
#' @export
`polarord` <-
    function(d, k=2)
{
    ## Create a zero matrix of ordination scores
    N <- attr(d, "Size")
    axes <- matrix(0, N, k)
    colnames(axes) <- paste0("PO", seq_len(k))
    rownames(axes) <- attr(d, "Labels")
    ## Get the total inertia: we divide with N to be consistent with
    ## wcmdscale, cmdscale, dbrda and and capscale
    inertia <- sum(d^2)/N
    ## return eigenvalues and endpoints for each axis
    ev <- numeric(k)
    names(ev) <- colnames(axes)
    endpoints <- matrix(0, 2, k,
                        dimnames = list(c("p1","p2"), paste0("PO", seq_len(k))))
    ## Iterate through dimensions 1..k
    for(dim in seq_len(k)) {
        m <- as.matrix(d)
        ## Find the first endpoint (variance method)
        p1 <- which.max(apply(m, 2, var))
        ## Find the second endpoint (regression method: regression
        ## coefficient was originally suggested, but covariance gives
        ## the same results, but faster)
        p2 <- which.min(cov(m[p1,],m))
        ## Project all points between these endpoints
        sco <- (m[p1,p2]^2  + m[p1,]^2 - m[p2,]^2)/2/m[p1,p2]
        ## Find the residual dissimilarities, with a guard against
        ## negative residuals in semimetric indices & oblique axes
        d <- sqrt(pmax(d^2 - dist(sco)^2,0))
        ## save the eigenvalue
        ev[dim] <- sum(dist(sco)^2)/N
        axes[,dim] <- sco
        endpoints[,dim] <- c(p1,p2)
    }
    out <- list(points = axes, inertia = inertia, eig = ev,
                endpoints = endpoints, call = match.call())
    class(out) <- "polarord"
    out
}

### polarord methods

#' @export
`print.polarord` <-
    function(x, ...)
{
    cat("Polar Ordination\n")
    cat("Call:", deparse(x$call), "\n\n")
    cat("Axis endpoints:\n")
    print(x$endpoints)
    cat("\nEigenvalues:\n")
    print(x$eig)
    cat("Total inertia:", x$inertia, "\n\n")
    if (!is.null(x$species))
        cat(gettextf("Species scores are expanded weighted averages from '%s'\n\n",
            attr(x$species, "data")))
}

#' @importFrom vegan ordiplot
#' @param x \code{polarord} result.
#' @param choices Axes shown.
#' @param type Type of graph which may be \code{"t"} for text, \code{"p"}
#'   for points or \code{"n"} for none (an empty plot).
#' @param display Items displayed: \code{"sites"} are always available,
#' but \code{"species"} only if they were added with sppscores.
#' @param \dots Other arguments to the function (passed to
#'   \code{\link[vegan]{ordiplot}}).
#' @rdname polarord
#'
#' @export
`plot.polarord` <-
    function(x, choices = c(1, 2), type = "t", display, ...)
{
    if (missing(display))
        if (is.null(x$species))
            display <- "sites"
        else
            display <- c("sites", "species")
    ordiplot(x, display = display, choices = choices, type = type, ...)
}

#' @importFrom vegan "sppscores<-" wascores
#' @export "sppscores<-"
#' @aliases sppscores<-
#' @param object \code{polarord} result.
#' @param value Community data to find the species scores.
#' @details Function \code{sppscores} can be used to add species scores
#'   to the ordination result.
#' @rdname polarord
#' @export
`sppscores<-.polarord` <-
        function(object, value)
{
    wa <- wascores(object$points, value, expand = TRUE)
    attr(wa, "data") <- deparse(substitute(value))
    object$species <- wa
    object
}

#' @importFrom vegan eigenvals
#' @export eigenvals
#' @aliases eigenvals
#' @rdname polarord
#' @export
`eigenvals.polarord` <-
    function(x, ...)
{
    out <- x$eig
    attr(out, "sumev") <- x$inertia
    class(out) <- "eigenvals"
    out
}
