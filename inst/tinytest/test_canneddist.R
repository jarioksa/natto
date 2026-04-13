require(vegan, quietly=TRUE)
data(spurn)
## example of one distance
expect_equivalent(canneddist(spurn, "steinhaus"),
                  vegdist(spurn, "bray"))
## check syntax of all canned distances
idx <- canneddist(help=TRUE)
expect_silent(for(nm in rownames(idx)) canneddist(spurn, nm))
