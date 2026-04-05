### distMeans gives distances of all points to the mean of xy columns
xy <- matrix(c(3,7,1,0,5,0,10,3,8,2), 5, 2)
expect_equivalent(as.matrix(dist(rbind(xy, colMeans(xy))))[6,-6],
                  distMeans(dist(xy)))
### distMeans will be at the centroid of PCoA
data(spurn)
d <- canneddist(spurn, "whittaker")
m <- cmdscale(distMeans(d, addcentre=TRUE), k = 2)
expect_equal(m[1,], c(0,0))
