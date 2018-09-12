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
        ## (A+B-2*J)/(A+B) with various names (and people rave on these)
        "sorensen" = list(method = "(A+B-2*J)/(A+B)", terms = "binary"),
        "bray" =   list(method = "(A+B-2*J)/(A+B)", terms = "minimum"),
        "whittaker" =  list(method = "(A+B-2*J)/(A+B)", terms = "binary"),
        "steinhaus" = list(method = "(A+B-2*J)/(A+B)", terms = "minimum"),
        ## Another popular one
        "jaccard" = list(method = "(A+B-2*J)/(A+B-J)", terms = "binary"),
        "ruzicka" = list(method = "(A+B-2*J)/(A+B-J)", terms = "minimum"),
        "similarityratio" = list(method = "(A+B-2*J)/(A+B-J)", terms = "quadratic"),
        ## Legendre & Legendre: metric distances
        "euclidean" = list(method = "sqrt(A+B-2*J)", terms = "quadratic"),
        "chord" = list(method = "sqrt(2*(1-J/sqrt(A*B)))", terms = "quadratic"),
        "geodesic" = list(method = "acos(J/sqrt(A*B))", terms = "quadratic"),
        "manhattan" = list(method = "A+B-2*J", terms = "minimum"),
        "ochiai" = list(method = "1-J/sqrt(A*B)", terms = "binary"),
        "cosine" = list(method = "1-J/sqrt(A*B)", terms = "quadratic"),
        ## Legendre & Legendre have some oddities, here as dissimilarities
        "triplejaccard" = list(method = "(A+B-2*J)/(A+B+J)", terms = "binary"),
        "sokalsneath" = list(method = "2*(A+B-2*J)/(2*A+2*B-3*J)", terms="binary"),
        "russellrao" = list(method = "1-J/P", terms = "binary"),
        "alt.kulczynski" = list(method = "(A+B-3*J)/(A+B-2*J)", terms = "binary"),
        ## Kulczynskis
        "kulczynski" = list(method =  "1-J/2/A-J/2/B", terms = "minimum"),
        "bin.kulczynski" = list(method = "1-J/2/A-J/2/B", terms="binary"),
        ## Raup-Crick with equal sampling probs
        "raup" = list(method = "1-phyper(J-1,A,P-A,B)", terms="binary"),
        ## 1 if no shared species, 0 if there is a shared species
        "shared" = list(method = "J==0", terms = "binary"),
        ## Hubalek 1982, Biol Rev 57, 669-689 adds (as distances):
        "braunblanquet" = list(method = "1-J/pmax(A,B)", terms="binary"),
        "simpson" = list(method = "1-J/pmin(A,B)", terms="binary"),
        "sorgenfrei" = list(method = "1-J*J/A/B", terms="binary"),
        "fager.mcgowan" = list(method = "0.5-J/sqrt(A*B)+pmax(A,B)/2",
                               terms="binary"),
        "mountford.init" = list(method = "2*J/(2*A*B-(A+B)*J)", terms="binary")
        )

    ind <- match.arg(method, names(index))
    z <- index[[ind]]
    designdist(x, method = z$method, terms = z$terms, name = ind)
}
