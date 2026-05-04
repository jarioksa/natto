### gnlmod: Arrhenius model of species accumulation.
require(vegan, quietly=TRUE) # data, SSarrhenius
data(sipoo, sipoo.map, package = "vegan")
S <- specnumber(sipoo)

expect_silent(mod <- gnlmod(S ~ SSarrhenius(area, k, z),
                            data=cbind("S" = S, sipoo.map),
                            family=poisson))
mod0 <- glm(S ~ log(area), sipoo.map, family=poisson)

expect_equal(deviance(mod), deviance(mod0))
expect_equivalent(fitted(mod), fitted(mod0))
k <- coef(mod0)
expect_equivalent(coef(mod), c(exp(k[1]), k[2]))
