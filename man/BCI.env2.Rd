\name{BCI.env2}
\alias{BCI.env2}
\docType{data}
\title{
  Habitat Types of Barro Colorado Island Forest Plots
}
\description{

  Habitat types of Barro Colorado Island Forest Dynamics Plots
  classified in 25 grid cells in each 1-ha plot (Harms et al. 2001).

}

\usage{data("BCI.env2")}

\format{

  A data frame with 50 observations on the area (in ha) of following
  habitat types:
  
  \describe{
    \item{\code{Young}}{Young forests of ca. 100 years.}
    \item{\code{Stream}}{Streamside vegetation.}
    \item{\code{Swamp}}{Swamp.}
    \item{\code{OldHigh}}{Old forests on flat surface above 152 m altitude.}
    \item{\code{OldLow}}{Old forests on flat surface at or below 152 m
       altitude.}
    \item{\code{OldSlope}}{Old forests on steep slopes (> 7 degrees).}
    \item{\code{Mixed}}{Mixed habitat.}
  }
}

\details{

  The data set provides habitat type classification for the Barro
  Colorado Island data \code{\link[vegan]{BCI}} in the \pkg{vegan}
  package. The \pkg{vegan} package provides alternative environmental
  data \code{\link[vegan]{BCI.env}} which has variables derived from
  this data set: Dominant habitat type excluding \code{Stream},
  presence of \code{Stream} grid cells in the plot, and habitat type
  diversity of each plot (see Examples).
  
  The data were extracted from Harms et al. (2001: Fig. 1) who give a
  map of 50 forest plots divided into 25 grid cells, and the areas are
  expressed in ha (at 0.04 ha resolution). 

}

\source{
  Harms et al. (2001), Fig. 1.
}

\references{

  Harms K.E., Condit R., Hubbell S.P. & Foster R.B. (2001) Habitat
  associations of trees and shrubs in a 50-ha neotropical forest
  plot. \emph{J. Ecol.} 89, 947--959.

}
\examples{

## Derive variables in BCI.env data in vegan
data(BCI.env2)
## 1. The dominant habitat type in each plot
k <- apply(BCI.env2[,-2], 1, which.max)
Habitat <- (names(BCI.env2)[-2])[k]
## 2. Presence of streamside habitats in the plot
Stream <- with(BCI.env2, ifelse(Stream > 0, "Yes", "No"))
## 3. Habitat type diversity (needs vegan)
if (require(vegan)) {
EnvHet <- diversity(BCI.env2, "simpson")
}

## The combined variables in BCI.env and the original area variables in
##  BCI.env2 give nearly equal constraints.
if (require(vegan)) {
data(BCI)
## habitat areas per 1-ha plot
ord <- cca(BCI ~ ., BCI.env2)
## dominant habitat type of each 1-ha plot
ord0 <- cca(BCI ~ Habitat)
plot(procrustes(ord, ord0))
}

} \keyword{datasets}
