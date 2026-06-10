#' Cast vegan::rda Result to stats::prcomp
#'
#' Function casts a result object of unconstrained
#' \code{\link[vegan]{rda}} or \code{\link[vegan]{pca}} to a
#' \code{\link{prcomp}} result object, but without the last zero
#' eigenvalue and rubbish eigenvectors of \code{prcomp}. The
#' translation is best evaluated by looking at the code of the
#' function.
#'
#' @param x An unconstrained \code{\link[vegan]{rda}} result object.
#'
#' @return A \code{\link{prcomp}} result object.
#'
#' @author Jari Oksanen

#' @rdname as.prcomp
#' @export
`as.prcomp` <-
    function(x)
{
    UseMethod("as.prcomp")
}

#' @importFrom vegan scores
#' @rdname as.prcomp
#' @export
`as.prcomp.rda` <-
    function(x)
{
    if (!is.null(x$CCA) || !is.null(x$pCCA))
        stop("works only with unconstrained rda")
    structure(
        list(sdev = sqrt(x$CA$eig),
             rotation = x$CA$v,
             center = attr(x$Ybar, "scaled:center"),
             scale = if(!is.null(scl <- attr(x$Ybar, "scaled:scale")))
                         scl
                     else
                         FALSE,
             x = scores(x, display = "sites", scaling = 1,
                        choices = seq_len(x$CA$rank),
             const = sqrt(x$tot.chi * (nrow(x$CA$u)-1)))),
        class = "prcomp")
}
