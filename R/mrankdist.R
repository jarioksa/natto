#' Mean Rank Shift Between Pairs of Sites
#'
#' Function calculates the mean rank shift between decreasing rank
#' orders of species (columns). The averaging is done by the number of
#' species occurring in both compared sites, and the results are
#' returned as dissimilarities.
#'
#' @param x The input data.
#'
#' @param missing If \code{TRUE}, missing species (zero abundance)
#' will not be used in comparing rank shifts, but rank shifts are only
#' compared for species that occur in both compared SUs. If
#' \code{FALSE}, the missing species are given tied rank as last
#' species.
#'
#' @references Collins, S.L, Suding, K.N., Cleland, E.E., Batty, M.,
#' Pennings, S.C., Gross, K.L., Grace, J.B., Gough, L., Fargione,
#' J.E. & Clark, C.M. (2008). Rank clocks and plant community
#' dynamics. \emph{Ecology} 89, 3534--3541.
#'
#' @importFrom stats dist
#' @export
`mrankdist` <-
    function(x, missing = TRUE)
{
    ## If missing=TRUE, rank shifts are only recorded for species
    ## present in both compared sites: make missing species NA
    if (missing) {
        x[x==0] <- NA
        x <- t(apply(-x, 1, rank, na.last = "keep"))
    } else {
        x <- t(apply(-x, 1, rank, na.last = TRUE))
    }
    d <- dist(x, method = "manhattan")
    ## Manhattan in basic dist() is divided by J/S where J is the
    ## number of species occurring in both sites (when missing species
    ## are NA), and S is the total number of species. We only want the
    ## average per J.
    d <- d/ncol(x)
    attr(d, "method") <- "meanrankshift"
    d
}
