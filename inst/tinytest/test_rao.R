### Tests for Rao diversity, distance & standardization

require(vegan, quietly=TRUE)
data(BCI, package = "vegan")
## Limiting cases: if all species are distinct, identical to basic methods.

## Simpson diversity:
d <- as.dist(matrix(1, ncol(BCI), ncol(BCI)))
expect_equivalent(qrao(BCI, d),
                  diversity(BCI, "simpson"))
## Euclidean distances
expect_equivalent(distrao(BCI, d, "euclidean"),
                  dist(decostand(BCI, "tot")))

## Use taxonomic distance
data(BCI.taxon)
d <- taxa2dist(BCI.taxon, varstep = TRUE)

## Rao standardization
expect_silent(Z <- raostand(BCI, d))
expect_equivalent(1 - rowSums(decostand(BCI, "tot")^2),
                  diversity(BCI, "simpson")) # simple Simpson
expect_equivalent(1 - rowSums(Z^2),
                  qrao(BCI, d)) # Rao quadratic entropy
expect_equivalent(1 - diag(tcrossprod(Z)),
                  qrao(BCI, d)) # Rao quadratic entropy
expect_equivalent(dist(Z)^2/2,
                  distrao(BCI, d)) # Rao's Jensen distance

