### taxondist: numerical example from Primer manual posted by user
### ailich in vegan issue #430.

require(vegan, quietly=TRUE) # taxa2dist
## Two sampling units, ten species
dat <- as.data.frame(rbind(c(1,1,1,1,0,0,0,1,0,1),
                           c(1,0,0,1,0,1,0,0,0,1)))
names(dat)<- paste0("Spec", 1:10)

## Taxonomy data.frame
dat.taxon<- structure(list(
    Genus = c("G1", "G1", "G2", "G3", "G4", "G4",  "G4", "G5", "G6", "G6"),
    Family = c("F1", "F1", "F1", "F1", "F2", "F2", "F2", "F3", "F3", "F3"),
    Order = c("O1", "O1", "O1", "O1", "O1", "O1", "O1", "O2", "O2", "O2"),
    Class = c("C1", "C1", "C1", "C1", "C1", "C1", "C1", "C1", "C1", "C1")),
    class = "data.frame",
    row.names = c("Spec1", "Spec2", "Spec3", "Spec4", "Spec5", "Spec6",
                  "Spec7", "Spec8", "Spec9", "Spec10"))

## Taxonomic dissimilarities
taxdis<- as.matrix(taxa2dist(dat.taxon, varstep=FALSE))

## Reference values were given in three digits and multiplied by 100
## (but we use max 1 for dissimilarities).
expect_equivalent(as.numeric(taxondist(dat, taxdis, "gamma")),
                  0.2, tol=0.0005)
expect_equivalent(as.numeric(taxondist(dat, taxdis, "theta")),
                  0.198, tol=0.0005)
