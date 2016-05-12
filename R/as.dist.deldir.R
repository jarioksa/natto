#' Cast deldir Object to Distances Between Neighbours
#'
#' Function casts a Delaunay triangulation from
#' \code{\link[deldir]{deldir}} (\pkg{deldir} package) to a
#' \code{\link{dist}} object, with numeric values for distances
#' between neighbours and \code{NA} for other distances.
#'
#' @param m A \code{\link[deldir]{deldir}} result.
#' @param diag logical value indicating whether the diagonal of the distance
#'          matrix should be printed by \code{\link{print.dist}}.
#' @param upper logical value indicating whether the upper triangle of the
#'          distance matrix should be printed by \code{\link{print.dist}}.
#' @importFrom stats as.dist
#'
#' @export
`as.dist.deldir` <-
    function(m, diag = FALSE, upper = FALSE)
{
    n <- m$n.data
    mat <- matrix(NA, n, n)
    for(i in seq_len(nrow(m$delsgs))) {
        d <- (m$delsgs$x1[i] - m$delsgs$x2[i])^2 + (m$delsgs$y1[i] - m$delsgs$y2[i])^2
        j <- max(m$delsgs$ind1[i], m$delsgs$ind2[i])
        k <- min(m$delsgs$ind1[i], m$delsgs$ind2[i])
        mat[j,k] <- sqrt(d)
    }
    mat <- as.dist(mat, diag = diag, upper = upper)
    attr(mat, "call") <- match.call()
    attr(mat, "method") <- "deldir"
    mat
}
