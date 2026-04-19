### rdecorana tests. The function uses poor iterative method (power
### method) with tol=sqrt(.Machine$double.eps) or appoximately
### 1.49e-8, and it makes no sense to use stricter tol. Different
### starting values of iteration can give results that expect_equal
### considers different and errors accumulate toward later
### axes. Therefore we need to reduce tol in the following tests.

## results should be similar as with vegan::decorana.
## However, there are numerical differences and usually
## tolerance must be relaxed to pass tests.
require(vegan, quietly=TRUE)
data(dune)
## 1. default call implying rescaling with iresc=4
m0 <- vegan::decorana(dune)
expect_silent(m <- rdecorana(dune))
expect_equivalent(m0$evals.decorana, m$evals, tol=1e-7)
expect_equivalent(diag(cor(m0$rproj, m$rproj)), c(1,1,1,1), tol=1e-5)
expect_equivalent(diag(cor(m0$cproj, m$cproj)), c(1,1,1,1), tol=1e-5)

## 2. no rescaling, only detrending
m0 <- vegan::decorana(dune, iresc = 0)
expect_silent(m <- rdecorana(dune, iresc = 0))
expect_equivalent(m0$evals.decorana, m$evals, tol=1e-7)
expect_equivalent(diag(cor(m0$rproj, m$rproj)), c(1,1,1,1), tol=1e-5)
expect_equivalent(diag(cor(m0$cproj, m$cproj)), c(1,1,1,1), tol=1e-5)

## 3. Orthogonal CA
m0 <- vegan::decorana(dune, ira=1)
expect_silent(m <- rdecorana(dune, ira=1))
expect_equivalent(m0$evals, m$evals, tol=1e-8)
expect_equivalent(diag(cor(m0$rproj, m$rproj)), c(1,1,1,1), tol=1e-6)
expect_equivalent(diag(cor(m0$cproj, m$cproj)), c(1,1,1,1), tol=1e-6)

## gradrescale: should be equal to decorana(t(x), iresc=1) column score 1
## works also as a test for rescaling stage in decorana
m0 <- vegan::decorana(dune, iresc = 0)
u <- m0$rproj[,1]
expect_silent(uscaled <- gradrescale(u, dune))
## one rescaling of transposed data to get rescaled "species" score
mtrans <- vegan::decorana(t(dune), iresc = 1)
expect_equal(cor(mtrans$cproj[,1], uscaled), 1, tol=1e-5)
expect_equal(order(u), order(uscaled), tol=1e-5,
             info="rescaling is monotonous")
expect_equal(range(mtrans$rproj[,1]), range(uscaled), tol=1e-5,
             info="gradient length matches transposed rproj")
