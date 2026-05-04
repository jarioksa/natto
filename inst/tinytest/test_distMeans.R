### distMeans gives distances of all points to the mean of xy columns
xy <- matrix(c(3,7,1,0,5,0,10,3,8,2), 5, 2)
expect_equivalent(distMeans(dist(xy)),
                  as.matrix(dist(rbind(xy, colMeans(xy))))[6,-6])
### distMeans will be at the centroid of PCoA
data(spurn)
d <- canneddist(spurn, "whittaker")
m <- cmdscale(distMeans(d, addcentre=TRUE), k = 2)
expect_equal(m[1,], c(0,0))
### verification of distMeans.formula: subtract factor level means
### from every species after which rowSums of squared centred data
### gives the squared Euclidean distances to the group centroid.
data(dune, dune.env, package = "vegan")
fmeans <- tapply(dune, dune.env$Management, colMeans)
fcentred <- dune - do.call(rbind, fmeans)[dune.env$Management, ]
expect_equal(distMeans(dist(dune) ~ Management, dune.env),
             sqrt(rowSums(fcentred^2)))
## formula with constant gives the same distances as the default
## method
d <- canneddist(dune, "kulczynski")
expect_equal(distMeans(d ~ 1),
             distMeans(d))
