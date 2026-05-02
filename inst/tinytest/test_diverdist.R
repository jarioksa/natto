### diverdist & diverclust
require(vegan, quietly=TRUE)
data(spurn, package="natto")
## renyi=, hill=TRUE defines increase in number of species
expect_equivalent(diverdist(spurn, renyi=0, hill=TRUE),
                  designdist(spurn, "(A+B-2*J)/2", "binary"))
## no real way to verify that the results are correct
expect_silent(cl <- diverclust(spurn, renyi=Inf, equalize="renyi"))
expect_silent(cl0 <- diverclust(decostand(spurn, "max", MARGIN=1),
                                renyi=Inf, equalize="no"))
expect_inherits(cl, "hclust")
expect_equal(cl$merge, cl0$merge)
expect_equal(cl$height, cl0$height)
