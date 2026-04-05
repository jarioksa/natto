require(vegan, quietly=TRUE)
data(spurn)
expect_equivalent(canneddist(spurn, "steinhaus"),
                  vegdist(spurn, "bray"))
