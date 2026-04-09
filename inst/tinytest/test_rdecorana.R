### rdecorana tests

## results should be similar as with vegan::decorana.
## However, there are numerical differences and usually
## tolerance must be relaxed to pass tests.
require(vegan, quietly=TRUE)
data(dune)
## 1. default call implying rescaling with iresc=4
m0 <- vegan::decorana(dune)
m <- rdecorana(dune)
expect_equivalent(m0$evals.decorana, m$evals, tol=1e-7)
expect_equivalent(abs(diag(cor(m0$rproj, m$rproj))), c(1,1,1,1), tol=1e-5)
expect_equivalent(abs(diag(cor(m0$cproj, m$cproj))), c(1,1,1,1), tol=1e-5)

## 2. no rescaling, only detrending
m0 <- vegan::decorana(dune, iresc = 0)
m <- rdecorana(dune, iresc = 0)
expect_equivalent(m0$evals.decorana, m$evals, tol=1e-7)
expect_equivalent(abs(diag(cor(m0$rproj, m$rproj))), c(1,1,1,1), tol=1e-5)
expect_equivalent(abs(diag(cor(m0$cproj, m$cproj))), c(1,1,1,1), tol=1e-5)

## 3. Orthogonal CA
m0 <- vegan::decorana(dune, ira=1)
m <- rdecorana(dune, ira=1)
expect_equivalent(m0$evals, m$evals, tol=1e-8)
expect_equivalent(abs(diag(cor(m0$rproj, m$rproj))), c(1,1,1,1), tol=1e-7)
expect_equivalent(abs(diag(cor(m0$cproj, m$cproj))), c(1,1,1,1), tol=1e-7)
