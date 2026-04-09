### PCO of distconstrain should give same solutions as
### constrained/unconstrained components of dbrda
require(vegan, quietly=TRUE)
data(dune, dune.env)
d <- canneddist(dune, "chord")
mod <- dbrda(d ~ Management + Moisture, dune.env)
expect_silent(dcon <- distconstrain(d ~ Management + Moisture,
                      dune.env))
expect_equivalent(eigenvals(pco(dcon)),
                  eigenvals(mod, "constrained"))
expect_silent(dres <- distconstrain(d ~ Management + Moisture,
                      dune.env, residuals = TRUE))
expect_equivalent(eigenvals(pco(dres)),
                  eigenvals(mod, "unconstrained"))
