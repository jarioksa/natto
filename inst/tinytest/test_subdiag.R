### subdiag
x <- matrix(0, 5, 5)
x <- row(x) + col(x)/10
expect_equal(subdiag(x),
             c(2.1, 3.2, 4.3, 5.4))
expect_equal(subdiag(x),
             subdiag(as.dist(x)))
