### Tests for constrained NMDS: caxNMDS, cdisNMDS
data(mite, mite.env, package = "vegan")
d <- canneddist(mite, "geodesic")

### caxNMDS
expect_silent(caxNMDS(d ~ WatrCont + SubsDens + Shrub, mite.env))
## verification test: using non-constrained NMDS axes as constraints
## will reproduce NMDS. This would fail in cdisNMDS.
m0 <- vegan::metaMDS(d, trace = 0)
NAX <- as.data.frame(m0$points)
expect_silent(m <- caxNMDS(d ~ ., NAX))
expect_equal(vegan::procrustes(m0, m)$ss, 0)


### cdisNMDS: no proper tests, only see that it runs.
expect_silent(
    cdisNMDS(d ~ WatrCont + SubsDens + Shrub, mite.env))

### verification: limiting case with over-defined constraints equals NMDS
data(spurn)
dummy <- factor(seq_len(nrow(spurn))) # rank of constraints = no. of points
d <- canneddist(spurn, "bray")
u <- cmdscale(d, k = 2) # to guarantee same starting configuration
m0 <- vegan::monoMDS(d, u, k = 2) # non-constrained NMDS from cmdscale
mdis <- cdisNMDS(d ~ dummy)
max <- caxNMDS(d ~ dummy)
## mdis$sites were added by MDSaddpoints and have symmetric SS 2e-6
expect_equal(vegan::procrustes(m0, mdis$constraints, symmetric = TRUE)$ss, 0)
expect_equal(vegan::procrustes(m0, max)$ss, 0)

## check that optimization does not stop at starting values with
## duplicated constraints.
data(dune, dune.env, package = "vegan")
d <- canneddist(dune, "barkman")
## 12 distinct combinations of constraints for 20 points
expect_silent(m <- caxNMDS(d ~ Management + Moisture, dune.env))
expect_true(m$counts[1] > 1) # optim progressed
expect_identical(duplicated(scores(m)), duplicated(m$model.matrix))
