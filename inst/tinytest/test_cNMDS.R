### Tests for constrained NMDS: caxNMDS, cdisNMDS
data(mite, mite.env, package = "vegan")
d <- canneddist(mite, "geodesic")

### caxNMDS
expect_silent(caxNMDS(d ~ WatrCont + SubsDens + Shrub, mite.env))
## verification test: using non-constrained NMDS axes as constraints
## will reproduce NMDS.
m0 <- vegan::metaMDS(d, trace = 0)
NAX <- as.data.frame(m0$points)
expect_silent(m <- caxNMDS(d ~ ., NAX))
expect_equal(vegan::procrustes(m0, m)$ss, 0)


### cdisNMDS: no proper tests, only see that it runs.  NB! Directly
### using 'd' i formula fails as distconstrain cannot find 'd':
### embedding distconstrain in expect_silent() fails due to scoping
### issues
expect_silent(
    cdisNMDS(canneddist(mite, "chord") ~ WatrCont + SubsDens + Shrub,
             mite.env))
