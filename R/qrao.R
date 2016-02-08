#' Rao's Quadratic Entropy and Dissimilarity
#'
#' Rao's quadratic entropy (Rao 1982) is a generalization of Simpson
#' (or Gini-Simpson) diversity index to a situation where species are
#' not completely independent. Simpson's diversity index is the
#' probability that two random species are different, but in Rao's
#' index the magnitude of difference (\eqn{0..1}) is used as a weight
#' to these probabilities.
#'
#' @param x Community data. The diversity is estimated for the rows of
#' the matrix.
#' @param d Dissimilarities among species (columns). If some some
#' dissimilarities are \eqn{>1}, these are divided with the maximum
#' dissimilarity. If all \eqn{d = 1}, the result will be equal to
#' Simpson's index.
#'
#' @return Vector of Rao's quadratic entropy values.
#'
#' @references Rao, C.R. (1982) Diversity and dissimilarity
#' coefficients: a unified approach. \emph{Theoretical Population
#' Biology} 21, 24--43.
#'
#' @seealso There are other implementations of this function in
#' \R. Most notably function \code{\link[ade4]{divc}} in \pkg{ade4}.
#' 
#' @examples
#' if (require(vegan)) {
#' data(dune, dune.phylodis)
#' qrao(dune, dune.phylodis)
#' }
#' 
#' @importFrom vegan decostand
#'
#' @export
`qrao` <-
    function(x, d)
{
    ## handle community matrix
    x <- as.matrix(x)
    x <- decostand(x, "tot")
    ## dissimilarities among species columns. These should be in 0..1,
    ## where 1 imeans completely different species, and 0 completely
    ## identical.
    d <- as.dist(d)
    if (max(d) > 1)
        d <- d/max(d)
    ## diversities
    n <- nrow(x)
    p <- ncol(x)
    ltri <- lower.tri(matrix(0, p, p))
    div <- numeric(n)
    for(i in seq_len(n))
        div[i] <- 2*sum(outer(x[i,], x[i,])[ltri] * d)
    div
}

#' @rdname qrao
#' @export
`distrao` <-
    function(x, d)
{
    x <- as.matrix(x)
    d <- as.matrix(as.dist(d))
    H <- x %*% d %*% t(x)
    diaH <- diag(H)
    out <- H - outer(diaH, diaH, "+")/2
    as.dist(out)
}
