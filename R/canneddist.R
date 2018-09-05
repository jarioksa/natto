#' Canned Dissimilarities with their Vernacular Names
#'
#' Function is a storehouse for dissimilarity indices that can be
#' called by their vernacular names. The function is based on
#' \code{\link[vegan]{designdist}} (\CRANpkg{vegan} package).
#'
#' @param x Input data
#' @param method Vernacular name for a dissimilarity index
#'
#' @return Function returns a \code{\link{dist}} object of dissimilarities.
#'
#' @author Jari Oksanen
#'
#' @seealso \code{\link[vegan]{designdist}}.
#'
#' @examples
#' data(spurn)
#' ## Ochiai dissimilarity
#' canneddist(spurn, "ochiai")
#'
#' @importFrom vegan designdist
#'
#' @export
`canneddist` <-
    function(x, method)
{
    index <- list(
        ## first proof-of-the-concept
        "sorensen" = list(method = "(A+B-2*J)/(A+B)", terms = "binary"),
        "bray" =   list(method = "(A+B-2*J)/(A+B)", terms = "minimum"),
        "whittaker" =  list(method = "(A+B-2*J)/(A+B)", terms = "binary"),
        "ochiai" = list(method = "1-J/sqrt(A*B)", terms = "binary"),
        "cosine" = list(method = "1-J/sqrt(A*B)", terms = "quadratic"),
        ## Legendre & Legendre: metric distances
        "euclidean" = list(method = "sqrt(A+B-2*J)", terms = "quadratic"),
        "chord" = list(method = "sqrt(2*(1-J/sqrt(A*B)))", terms = "quadratic"),
        "geodesic" = list(method = "acos(J/sqrt(A*B))", terms = "quadratic"),
        "manhattan" = list(method = "A+B-2*J", terms = "minimum"),
        ## L & L: binary similarity indices as distances 1-S in terms A,B,J
        "jaccard" = list(method = "(A+B-2*J)/(A+B-J)", terms = "binary"),
        "triplejaccard" = list(method = "(A+B-2*J)/(A+B+J)", terms = "binary"),
        "sokalsneath" = list(method = "2*(A+B-2*J)/(2*A+2*B-3*J)", terms="binary"),
        "russellrao" = list(method = "1-J/P", terms = "binary"),
        "bin.kulczynski14" = list(method = "(A+B-3*J)/(A+B-2*J)", terms = "binary"),
        "bin.kulczynski15" = list(method = "1-J/2/A-J/2/B", terms="binary"),
        ## L & L: quantitative similarity indices as distances 1-S
        "steinhaus" = list(method = "(A+B-2*J)/(A+B)", terms = "minimum"),
        "kulczynski" = list(method =  "1-J/2/A-J/2/B", terms = "minimum"))

    ind <- match.arg(method, names(index))
    z <- index[[ind]]
    designdist(x, method = z$method, terms = z$terms, name = ind)
}
