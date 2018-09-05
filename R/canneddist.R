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
        "sorensen" = list(method = "(A+B-2*J)/(A+B)", terms = "binary"),
        "bray" =   list(method = "(A+B-2*J)/(A+B)", terms = "minimum"),
        "whittaker" =  list(method = "(A+B-2*J)/(A+B)", terms = "binary"),
        "ochiai" = list(method = "1-J/sqrt(A*B)", terms = "binary"),
        "cosine" = list(method = "1-J/sqrt(A*B)", terms = "quadratic"))
    ind <- match.arg(method, names(index))
    z <- index[[ind]]
    designdist(x, method = z$method, terms = z$terms, name = ind)
}
