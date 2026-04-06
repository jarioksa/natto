## test for as.prcomp.rda structure
require(vegan, quietly=TRUE) # rda
data(spurn)

mrda <- pca(spurn) # or equivalently rda(spurn)

expect_equal(names(prcomp(spurn)), names(as.prcomp(mrda)))
expect_inherits(as.prcomp(mrda), "prcomp")
