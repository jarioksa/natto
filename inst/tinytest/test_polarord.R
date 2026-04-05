### polar ordination

## Recovers Euclidean coordinates
require(vegan, quietly = TRUE) # procrustes
set.seed(4711)
X <- as.data.frame(matrix(runif(30*3), nrow=30))
d <- dist(X) # Euclidean distances
m <- polarord(d, k=3)
## finds coordinates in X: Procrustes SS nearly zero
expect_equal(procrustes(polarord(d, k=3), X)$ss, 0)
