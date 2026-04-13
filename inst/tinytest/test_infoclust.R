### infoclust: replicate numerical example in Legendre & Legendre
### (2012) Numerical Ecology, p. 372-376.
data(pond)
expect_silent(cl <- infoclust(pond))
## pairwise information distances in three decimals (p. 374 top)
expect_equal(as.numeric(canneddist(pond, "information")),
             c(2.773, 8.318, 9.704, 9.704,
               8.318, 9.704, 6.931,
               4.159, 4.159,
               2.773), tol = 0.0005)
## Fusion levels given in three decimals (p. 374 bottom)
expect_equal(sort(unique(cl$height)),
             c(2.773, 7.638, 26.920), tol=0.0005)
## Same clusters (Fig. 8.15, p. 375)
## > cutree(cl, 3)
## 212 214 233 431 432
##   1   1   2   3   3
expect_equal(cutree(cl, 3),
             structure(c(1,1,2,3,3), names = rownames(pond)))
