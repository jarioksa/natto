## Welcome to natto!

[**natto**](http://en.wikipedia.org/wiki/Nattō) is an extreme
**vegan** package that contains experimental, specialized or
unfinished code that does not fit to the more general **vegan**
package. Some of these functions may be merged to **vegan** later, but
the package is not a code kindergarten with guaranteed tenure track to
**vegan**.

Currently **natto** contains

 * Hierarchical clustering based on information analysis after Williams,
   Lambert & Lance, *J. Ecol.* **54,** 427-445 (1966).
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
