#' Return Gammaplus and Thetaplus Elements of Clarke's Taxonomic Community Dissimilarity
#'
#' @param x Community data; will be treated as binary presence/absence matrix.
#' @param taxdist Taxonomic, phylogenetic or other dissimilarity matrix among species (columns of x).
#'
#' @return list of Clarke's Gammaplus and Thetaplus elements as distance structures.
vegtaxdist <- 
    function (x, taxdist) 
{
    x <- as.matrix(x)
    x <- ifelse(x > 0, 1, 0)
    taxdist <- as.matrix(taxdist)
    if (NCOL(taxdist) != NCOL(x)) stop("Number of columns do not match in 'x' and 'taxdist'")
    N <- NROW(x)
    Gammaplus <- matrix(0, N, N)
    Thetaplus <- matrix(0, N, N)
    for(j in 1:(N-1)) {
        for(i in (j+1):N) {
           crosstd <- (outer(x[i,], x[j,]) * taxdist)[x[i,] > 0, x[j,] > 0, drop = FALSE]
           min1 <- apply(crosstd, 1, min)
           min2 <- apply(crosstd, 2, min)
           Gammaplus[i,j] <- mean(c(min1, min2))
           Thetaplus[i,j] <- (mean(min1) + mean(min2))/2
        }
    }
    list(Gammaplus = as.dist(Gammaplus), Thetaplus = as.dist(Thetaplus))
}
