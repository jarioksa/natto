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
