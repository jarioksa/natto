### Position and Species Vector Ordinations

## verification: replicate numerical example in Orloci, Handbook of
## Vegetation Science 5, 249-286 (1973).

## Orloci's input matrix X is transposed to vegan standard: 4 columns
## for quadrats, 5 rows for species (p. 267).
X <- matrix(c(1,3,4,4,
              3,3,1,1,
              5,3,5,3,
              1,1,3,3,
              4,2,1,1),
            nrow=5, ncol=4, byrow = TRUE)
## Numerical results are based on cross-products instead of variances
## natto uses: we need to adjust by degrees of freedom n-1 (or 4-1=3)
## for eigenvalues, and its sqrt for eigenvectors. Results are given
## in 3 decimals on p. 272-273. First eigenvalues
S <- c(17.818, 5.182, 1.000)
## Eigenvectors: first is exact, others in 3 decimals, and the last
## one is reversed to natto.
pvorloci <- cbind("PVO1" = c(11, 1, -5, -7)/sqrt(11),
                  "PVO2" = c(0, 1.706, -1.492, -0.213),
                  "PVO3" = -c(0, 0, -0.707, 0.707))
## tests
expect_silent(mod <- posvectord(t(X))) # transposed
expect_equivalent(mod$eig * 3, S, tol=0.0005) # 3 decimals
expect_equal(mod$points[,1] * sqrt(3), pvorloci[,1]) # exact
expect_equal(mod$points * sqrt(3), pvorloci, tol=0.0005) # 3 decimals

## verification: recover Euclidean coordinates
require(vegan, quietly = TRUE) # procrustes
set.seed(4711)
X <- as.data.frame(matrix(runif(40*4), nrow=40))
## posvectord uses Euclidean algebra and recover original coordinates even
## when the axis are forced to go through sample points
expect_equal(procrustes(posvectord(X)$points, X)$ss, 0)

## spvectord: Orlóci 1978 Multivariate analysis in vegetation research
## (2 ed.), p. 28-31 gives numerical example without specifying the
## name of the method. Some results are given in 7 digits, but it
## seems that the calculations were made on single precision and only
## 5 to 6 digits are exact. The numerical equality in R uses tolerance
## based on machine epsilon, or in double precision eps = 2^-52 and
## tol = sqrt(eps) = 1.5e-8. In single precision eps = 2^-23 and
## tol = sqrt(eps) = 0.00035, and we need to use slacker tol in
## comparing Orlóci numbers to spvectord results.

## Table 1-1, page 13
tab11 <- matrix(c(45,18,3,10,9,
                  2,92,40,61,53,
                  9,32,5,11,99,
                  26,48,83,3,21,
                  3,73,68,32,49,
                  5,80,27,2,81,
                  2,95,2,39,72,
                  90,13,17,2,6,
                  31,92,1,8,90,
                  16,78,23,6,62),
                nrow=5, ncol=10, byrow=FALSE)
dimnames(tab11) <- list(paste0("spe", 1:5),
                        paste0("quad", 1:10)) # transposed
## Orlóci results (p. 28-31), sums of squares (SS) and ranks of species
SS <- c(1514.55, 17027.20, 9061.53, 2956.14, 6281.31)
rank <- c(5, 1, 2, 4, 3) # Species 2 is the most important
SS <- SS[order(rank)] # sort SS in decreasing order
## natto::spvectord
TOL <- sqrt(2^-23)/100
expect_silent(mod <- spvectord(t(tab11)))  # t(tab11)
expect_equivalent(mod$eig * 9, SS, tol = TOL) # natto uses var, Orlóci SS
expect_equivalent(mod$totvar * 9, sum(SS), tol = TOL)

## NB, Orlóci 1973 _Nature_ paper has a numerical example, but the
## data are incorrect for variable 'Style': its mean and variance do
## not match the observed numbers.
