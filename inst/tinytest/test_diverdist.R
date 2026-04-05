### diverdist & diverclust
require(vegan, quietly=TRUE)
data(spurn)
## renyi=, hill=TRUE defines increase in number of species
expect_equivalent(diverdist(spurn, renyi=0, hill=TRUE),
                  designdist(spurn, "(A+B-2*J)/2", "binary"))
## no real way to verify that the results are correct
cl <- diverclust(spurn)
expect_inherits(cl, "hclust")
