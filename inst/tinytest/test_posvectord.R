### Position and Species Vector Ordinations

## Recover Euclidean coordinates
require(vegan, quietly = TRUE) # procrustes
set.seed(4711)
X <- as.data.frame(matrix(runif(40*4), nrow=40))
## posvectord uses Euclidean algebra and recover original coordinates even
## when the axis are forced to go through sample points
expect_equal(procrustes(posvectord(X), X)$ss, 0)
## no verification test for spvectord
