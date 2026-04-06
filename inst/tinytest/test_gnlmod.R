### gnlmod: Arrhenius model of species accumulation.
require(vegan, quietly=TRUE) # data, SSarrhenius
data(sipoo, sipoo.map)
S <- specnumber(sipoo)

mod <- gnlmod(S ~ SSarrhenius(area, k, z), data=cbind("S" = S, sipoo.map),
              family=poisson)
mod0 <- glm(S ~ log(area), sipoo.map, family=poisson)

expect_equal(deviance(mod0), deviance(mod))
expect_equivalent(fitted(mod0), fitted(mod))
k <- coef(mod0)
expect_equivalent(c(exp(k[1]), k[2]), coef(mod))
