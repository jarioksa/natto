#' Compare Slopes of Mantel Regression Within and Across Blocks
#'
#' Functions fit linear regression on distances within blocks and
#' across blocks and compere their slopes using permutation tests.
#'
#' Functions are intended to analyse the differences of linear Mantel
#' regression within and across blocks. Typically the dependent
#' dissimilarities (\code{ydist}) are based on community data,
#' independent dissimilarities (\code{xdist}) are based on
#' environmental data, and \code{block} defines subsets of data, such
#' as geographical subareas.
#' 
#' Functions \code{slopediff} and \code{slopediff2} differ in their
#' permutation models. Function \code{slopediff} permutes the
#' classification vector defining \code{block}s, and keeps the
#' dependent (\code{ydist}) and independent (\code{xdist}) pairs of
#' dissimilarities fixed. Function \code{slopediff2} permutes
#' dependent data, and keeps \code{block}s and dependent distances
#' fixed.
#' 
#' @author Jari Oksanen
#'

#' @param ydist,xdist Dependent and indpendent dissimilarities.
#' @param block Blocking of sites.
#' @param nperm Number of simple permutations.
#'
#' @importFrom stats as.dist coef lm
#' @rdname slopediff
#' @export
`slopediff` <-
    function (ydist, xdist, block, nperm = 999) 
{
    within <- as.dist(outer(block, block, "=="))
    m1 <- lm(ydist ~  xdist, subset = within==1)
    m2 <- lm(ydist ~ xdist, subset = within==0)
    out <- matrix(NA, nrow=nperm+1, 4)
    out[nperm+1, ] <- c(coef(m1), coef(m2))
    for(i in 1:nperm) {
       k <- sample(length(block))
       w <- as.dist(outer(block[k], block[k], "=="))
       m1 <- lm(ydist~ xdist, subset=w==1)
       m2 <- lm(ydist ~xdist, subset = w==0)
       out[i,] <- c(coef(m1), coef(m2))
   }
   class(out) <- "slopediff"
   out
}
#' @importFrom stats as.dist coef lm
#' @rdname slopediff
#' @export
`slopediff2` <-
    function (ydist, xdist, block, nperm = 999) 
{
    within <- as.dist(outer(block, block, "=="))
    m1 <- lm(ydist ~  xdist, subset = within==1)
    m2 <- lm(ydist ~ xdist, subset = within==0)
    out <- matrix(NA, nrow=nperm+1, 4)
    out[nperm+1, ] <- c(coef(m1), coef(m2))
    ymat <- as.matrix(ydist)
    for(i in 1:nperm) {
       k <- sample(length(block))
       y <- as.dist(ymat[k,k])
       m1 <- lm(y ~ xdist, subset=within==1)
       m2 <- lm(y ~ xdist, subset = within==0)
       out[i,] <- c(coef(m1), coef(m2))
   }
   class(out) <- "slopediff"
   out
}
#' @param x Result output from \code{slopediff} or \code{slopediff2}.
#' @param col Plotting colour.
#' @param obsline Draw line on observed statistic in the
#' \code{\link{density}} plot. .
#'
#' @importFrom graphics plot abline
#' @importFrom stats density
#' 
#' @rdname slopediff
#' @export
`plot.slopediff` <-
    function (x, col=1, obsline = TRUE, ...) 
{
    add <- FALSE
    op <- par(mfrow=c(2,1))
    on.exit(par(op))
    if (!add) {
    	plot(density(x[,1] - x[,3]), col = col, main="Intercept", ...)
      if(obsline) abline(v=x[nrow(x),1] - x[nrow(x),3])
      plot(density(x[,2] - x[,4]), col = col, main="Slope", ...)
      if(obsline) abline(v=x[nrow(x),2] - x[nrow(x),4])
    } 
}
#' @param object Result object from \code{slopediff} or \code{slopediff2}.
#'
#' @rdname slopediff
#' @export
`summary.slopediff` <-
    function (object, ...) 
{
   nperm <- nrow(object)-1
   cat("permutation test based on", nperm, "permutations\n")
   tab <- matrix(0, nrow=4, ncol=2)
   tab[1,] <- object[nperm+1,1:2]
   tab[2,] <- object[nperm+1,3:4]
   tab[3,] <- tab[1,] - tab[2,]
   sig <- abs(cbind(object[,1] - object[,3], object[,2]-object[,4]))
   tab[4,] <- (colSums(sweep(sig, 2, sig[nperm+1,], ">"))+1)/(nperm+1)  
   colnames(tab) <- c("Intercept", "Slope")
   rownames(tab) <- c("within", "across", "difference", "P")
   tab
}
