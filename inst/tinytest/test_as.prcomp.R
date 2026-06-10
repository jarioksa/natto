## test for as.prcomp.rda structure
require(vegan, quietly=TRUE) # rda, procrustes
data(spurn)

mrda <- pca(spurn) # or equivalently rda(spurn)
mprc <- prcomp(spurn)

expect_silent(mcast <- as.prcomp(mrda))
expect_inherits(as.prcomp(mrda), "prcomp")
expect_equal(names(mprc), names(mcast))
## prcomp returns one axis beyond rank with "zero" eigenvalue (in
## practice below 1e-8) and rubbish eigenvectors, and we need to drop
## that from comparison
rank <- seq_len(mrda$CA$rank)
expect_equal(mcast$center, mprc$center)
expect_equivalent(mcast$sdev, mprc$sdev[rank])
expect_equal(procrustes(mcast$rotation, mprc$rotation[,rank])$ss, 0)
expect_equal(procrustes(mcast$x, mprc$x[,rank])$ss, 0)
