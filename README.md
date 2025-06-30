## Welcome to natto!

[**natto**](http://en.wikipedia.org/wiki/Nattō) is an extreme
**vegan** package that contains experimental, specialized or
unfinished code that does not fit to the more general **vegan**
package. Some of these functions may be merged to **vegan** later, but
the package is not a code kindergarten with guaranteed tenure track to
**vegan**.

Currently **natto** contains

 * Rao's quadratic entropy, phylogenetic dissimilarities and
   phylogenetic data standardization, plus related Clarke's taxonomic
   dissimilarity.
 * Polar Ordination A.K.A. Bray-Curtis Ordination.
 * Orlóci's Position Vector Ordination and a related method using
   species as ordination axes.
 * Implementation of Decorana in **R**: easier to study the algorithm
   and develop the method than with Fortran code.
 * Estimate residual species associations in CCA and RDA.
 * Randomized clustering and ordination based on Bayesian Jaccard
   dissimilarity with applications. Dissimilarities are random numbers
   from Beta distribution and each evaluation yields a new set of
   dissimilarities.  For clustering, this can be used to assess the
   support of clusters, and for ordination the confidence areas of
   coordinates.
 * Hierarchical clustering based on information analysis after Williams,
   Lambert & Lance, *J. Ecol.* **54,** 427-445 (1966).
 * beta and gamma diversity clustering.
 * Plotting cluster dendrograms (`hclust` results) so that terminal
   leaves are fans with base proportional to weight (size): `ginkgogram`.
 * Humped no-interaction model for species richness _vs._ biomass
   (transferred from the **vegan** package).
 * IAP functions to find the indicator values for species and statistics
   for sites in Index of Atmospheric Purity. The main focus is in indicator
   values for species which are the number of companion species. These can
   be used to find species that indicate high or low biodiversity, and 
   permutation test is provided to assess these values.
 * Mean Rank Shift function of Collins et al., *Ecology* **89,**
   3534-3541 (2008).
 * Testing differences of slopes in Mantel style regression between
   dissimilarities.
 * Casting unconstrained `vegan::rda()` result object to `stats::prcomp` 
   object.
 * Casting Delaunay triangulation from `deldir::deldir` result object to
   `stats::dist` object.

### Installation

Install **natto** directly from `github`:
```R
install.packages("devtools") # if you have not yet installed devtools
library(devtools)
install_github("jarioksa/natto")
```
Binary package can be installed via R-Universe:
```r
install.packages("natto", repos = c("https://jarioksa.r-universe.dev",
"https://cloud.r-project.org"))
```
