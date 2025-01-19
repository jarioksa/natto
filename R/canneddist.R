#' Canned Dissimilarities with their Vernacular Names
#'
#' Function is a storehouse for dissimilarity indices that can be
#' called by their vernacular names. The function is based on
#' \code{\link[vegan]{designdist}} (\CRANpkg{vegan} package).
#'
#' @details
#'
#' Function wraps popular dissimilarity indices for
#' \code{\link[vegan]{designdist}} allowing these to be called by
#' their popular names. It can have synonymous names for one
#' dissimilarity index. Use argument \code{help=TRUE} to see the
#' current selection of indices and their definitions.
#'
#' The function uses the main notation of
#' \code{\link[vegan]{designdist}} where terms are based on sums and
#' paired minima or sum of squares and crossproducts of pairs of
#' sampling units. For vectors \code{x} and \code{y} the
#' \code{"quadratic"} terms are \code{J = sum(x*y)}, \code{A =
#' sum(x^2)}, \code{B = sum(y^2)}, and \code{"minimum"} terms are
#' \code{J = sum(pmin(x,y))}, \code{A = sum(x)} and \code{B = sum(y)},
#' and \code{"binary"} terms are either of these after transforming
#' data into binary form (number of shared species \code{J}, and
#' number of species for each row, \code{A, B}.). Number of columns
#' (species) is notated as \code{P}, and the number of sampling units
#' is \code{N}.
#'
#' There is a huge number of indices, and the current selection is far
#' from comprehensive (but it can easily expanded). See References for
#' sources. Many sources use different notation, but they were changed
#' to the notation described above. For instance, in popular (but
#' strange) 2x2 contingency table notation for binary data \code{a =
#' J}, \code{b = A-J}, \code{c = B-J}, \code{d = P-A-B+J}. Some of
#' formulae may be surprising, but they are mathematically equivalent
#' to traditional ones. I challenge you to inspect Euclidean distance,
#' and once you see how it is derived, try Chord distance.
#'
#' @encoding UTF-8
#' @references
#'
#' Clarke, K.R., Somerfield, P.J. & Chapman, M.G. (2006).
#' On resemblance measures for ecological studies,
#' including taxonomic dissimilarities and a zero-adjusted
#' Bray-Curtis coefficient for denuded assemblages.
#' \emph{J. Exp. Marine Biol. & Ecol.} 330, 55-80. (\code{bray0}).
#'
#' \enc{Hubálek}{Hubalek}, Z. (1982). Coefficients of association and
#' similarity, based on binary (presence-absence) data: an
#' evaluation. \emph{Biological Review} 57, 669--689.
#'
#' Legendre, P. & Legendre, L. (2012). \emph{Numerical Ecology.} 3rd
#' English Ed., Elsevier.
#'
#' Yue, J.C. & Clayton, M.K. (2005). A similarity measure based on
#' species proportions. \emph{Communications in Statistics Theory and
#' Methods} 34, 2123--2131. \doi{10.1080/STA-200066418}.
#'
#' @param x Input data.
#' @param method Vernacular name for a dissimilarity index.
#' @param help List available indices and their definitions instead of
#'     calculating dissimilarities.
#'
#' @return Function returns a \code{\link{dist}} object of dissimilarities.
#'
#' @author Jari Oksanen
#'
#' @seealso
#'
#' Function is a wrapper to
#' \code{\link[vegan]{designdist}}. \CRANpkg{vegan} function
#' \code{\link[vegan]{betadiver}} is a similar collection of canned
#' indices for beta diversity, and many of these are well-known
#' dissimilarity indices.
#'
#' @examples
#' data(spurn)
#' ## Ochiai dissimilarity
#' canneddist(spurn, "ochiai")
#'
#' @importFrom stats phyper
#' @importFrom vegan designdist
#'
#' @export
`canneddist` <-
    function(x, method, help = FALSE)
{
    index <- list(
        ## (A+B-2*J)/(A+B) with various names (and people rave on these)
        "sorensen" = list(method = "(A+B-2*J)/(A+B)", terms = "binary",
                          maxdist = 1),
        "bray" =   list(method = "(A+B-2*J)/(A+B)", terms = "minimum",
                        maxdist = 1),
        "whittaker" =  list(method = "(A+B-2*J)/(A+B)", terms = "binary",
                            maxdist = 1),
        "steinhaus" = list(method = "(A+B-2*J)/(A+B)", terms = "minimum",
                           maxdist = 1),
        ## Zero-adjusted Bray Curtis of Clarke & al. (2006) J Exp Marine
        ## Biol & Ecol 330:55-80
        "bray0" = list(method = "(A+B-2*J)/(A+B+2*min(x[x>0]))",
                                 terms = "minimum", maxdist = NA),
        ## Another popular one
        "jaccard" = list(method = "(A+B-2*J)/(A+B-J)", terms = "binary",
                         maxdist = 1),
        "ruzicka" = list(method = "(A+B-2*J)/(A+B-J)", terms = "minimum",
                         maxdist = 1),
        "similarityratio" = list(method = "(A+B-2*J)/(A+B-J)",
                                 terms = "quadratic", maxdist = 1),

        ## Yue & Clayton (2005) Commun Stat Theory Methods 23, 2123-2131
        "yueclayton" = list(method = "(A+B-2*J)/(A+B-J)", terms = "quadratic",
                            maxdist = 1),
        ## Legendre & Legendre: metric distances
        "euclidean" = list(method = "sqrt(A+B-2*J)", terms = "quadratic",
                           maxdist = NA),
        "chord" = list(method = "sqrt(2*(1-J/sqrt(A*B)))", terms = "quadratic",
                       maxdist = sqrt(2)),
        "geodesic" = list(method = "acos(J/sqrt(A*B))", terms = "quadratic",
                          maxdist = pi/2),
        "ochiai" = list(method = "1-J/sqrt(A*B)", terms = "binary",
                        maxdist = 1),
        "cosine" = list(method = "1-J/sqrt(A*B)", terms = "quadratic",
                        maxdist = 1),
        "manhattan" = list(method = "A+B-2*J", terms = "minimum",
                           maxdist = NA),
        "information" = list(method = "log(4)*(A+B-2*J)", terms = "binary",
                             maxdist = NA),
        ## Legendre & Legendre have some oddities, here as dissimilarities
        "triplejaccard" = list(method = "(A+B-2*J)/(A+B+J)", terms = "binary",
                               maxdist = 1),
        "sokalsneath" = list(method = "2*(A+B-2*J)/(2*A+2*B-3*J)",
                             terms="binary", maxdist = 1),
        "russellrao" = list(method = "1-J/P", terms = "binary", maxdist = 1),
        ## Kulczynskis
        "kulczynski" = list(method =  "1-J/2/A-J/2/B", terms = "minimum",
                            maxdist = 1),
        "b.kulczynski" = list(method = "1-J/2/A-J/2/B", terms="binary",
                              maxdist = 1),
        ## Raup-Crick with equal sampling probs
        "raup" = list(method = "1-phyper(J-1,A,P-A,B)", terms="binary",
                      maxdist = 1),
        ## 1 if no shared species, 0 if there is a shared species
        "shared" = list(method = "J==0", terms = "binary", maxdist = 1),
        ## Hubalek 1982, Biol Rev 57, 669-689 adds (as distances):
        "braunblanquet" = list(method = "1-J/pmax(A,B)", terms="binary",
                               maxdist = 1),
        "simpson" = list(method = "1-J/pmin(A,B)", terms="binary", maxdist = 1),
        "sorgenfrei" = list(method = "1-J*J/A/B", terms="binary", maxdist = 1),
        "mountford.init" = list(method = "pmax(1-2*J/(2*A*B-(A+B)*J),0)",
                                terms="binary", maxdist = 1)
        )
    if (help)
        return(t(sapply(index, data.frame, stringsAsFactors = TRUE)))
    ind <- match.arg(method, names(index))
    z <- index[[ind]]
    dis <- designdist(x, method = z$method, terms = z$terms, name = ind,
                      maxdist = z$maxdist)
    attr(dis, "call") <- match.call()
    dis
}
