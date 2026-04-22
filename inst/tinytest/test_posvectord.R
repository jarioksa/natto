### Position and Species Vector Ordinations

## verification: replicate numerical example in Orloci, Handbook of
## Vegetation Science 5, 249-286 (1973).

## Orloci's input matrix X is transposed to vegan standard: 4 columns
## for quadrats, 5 rows for species (p. 257).
X <- matrix(c(1,3,4,4,
              3,3,1,1,
              5,3,5,3,
              1,1,3,3,
              4,2,1,1),
            nrow=5, ncol=4, byrow = TRUE)
## Numerical results are based on cross-products instead of variances
## natto uses: we need to adjust by degrees of freedom n-1 (or 4-1=3)
## for eigenvalues, and its sqrt for eigenvectors. Results are given
## in 3 decimals on p. 262-263. First eigenvalues
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
expect_equal(procrustes(posvectord(X), X)$ss, 0)
## spvectord: only test that it runs
data(spurn)
expect_silent(spvectord(spurn))
