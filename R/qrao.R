#' Rao's Quadratic Entropy
#'
#' Rao's quadratic entropy is a generalization of Simpson's (or Gini's
#' and Simpson's) diversity index to a situation where species are not
#' completely independent. Simpson's diversity index is the
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
