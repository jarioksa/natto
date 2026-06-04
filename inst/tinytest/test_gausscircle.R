## Verification test for gausscircle:
## Generate error-free 2D Gaussian response, and estimate its
## generating species parameters with gausscircle.

nobs <- 30
nsp <- 6
p <- matrix(0, nsp, 4)
colnames(p) <- c("xopt", "yopt", "tol", "top")
p[,1:2] <- runif(2*nsp, -2, 2)
p[,3] <- rnorm(nsp, 1, 0.2)
p[,4] <- rexp(nsp, 1/5)
x <- matrix(runif(2 * nobs, -2, 2), nrow = nobs)
y <- matrix(NA, nrow = nobs, ncol = nsp)
for(i in 1:nsp)
    y[,i] <- p[i,4] * exp(-rowSums(sweep(x, 2, p[i,1:2])^2)/2/p[i,3]^2)
## Re-find input response paramaters p
expect_silent(res <- gausscircle(x, y))
expect_equal(res[, 1:4], p)
## gausellipse: three first setup parameters should be equal, and in
## gaussellipse xtol == ytol (except names) and rxy == 0
expect_silent(res <- gaussellipse(x, y))
expect_equivalent(res[, 1:3], p[, 1:3])
expect_equivalent(res[, "ytol"], res[, "xtol"])
expect_equal(res[, "rxy"], rep.int(0, nsp))
## test for unit-tolerance responses
for (i in 1:nsp)
    y[,i] <- p[i,4] * exp(-rowSums(sweep(x, 2, p[i,1:2])^2)/2)
expect_silent(res <- gausscircle(x, y, unit = TRUE))
expect_equal(res[, c(1,2,4)], p[, c(1,2,4)])
expect_equal(gausscircle(x, y, unit = TRUE), gausscircle(x, y))

