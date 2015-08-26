#' Cast vegan::rda Result to base::prcomp
#'
#' Function casts a result object of unconstrained
#' \code{\link[vegan]{rda}} to a \code{\link{prcomp}} result object.
#'
#' @param x An unconstrained \code{\link[vegan]{rda}} result object.
#'
#' @importFrom vegan scores
#' @export
`as.prcomp.rda` <-
    function(x)
{
    if (!is.null(x$CCA) || !is.null(x$pCCA))
        stop("works only with unconstrained rda")
    structure(
        list(sdev = sqrt(x$CA$eig),
             rotation = x$CA$v,
             center = attr(x$CA$Xbar, "scaled:center"),
             scale = if(!is.null(scl <- attr(x$CA$Xbar, "scaled:scale")))
                         scl
                     else
                         FALSE,
             x = scores(x, display = "sites", scaling = 1,
             choices = seq_len(x$CA$rank),
             const = sqrt(x$tot.chi * (nrow(x$CA$u)-1)))),
        class = "prcomp")
}
