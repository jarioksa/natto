\name{BCI.taxon}
\alias{BCI.taxon}
\docType{data}
\title{
   Botanical Classification of Barro Colorado Island Species
}
\description{
  Classification of \code{\link[vegan]{BCI}} species based on
  APG IV (2016) and Chase & Reveal (2009).
}
\usage{data("BCI.taxon")}

\format{
  A data frame with 225 observations on the following 7 variables.
  \describe{
    \item{\code{genus}}{Genus.}
    \item{\code{family}}{APG IV family.}
    \item{\code{order}}{APG IV order.}
    \item{\code{superorder}}{Superorder (Chase & Reveal, 2009).}
    \item{\code{unrank}}{Either \code{eudicotyledons},
      \code{Magnoliidae} or \code{monocotyledons} which correspond to
      traditional major clades, but have no formal status in the
      classification of the subclass of Angiosperms. }
   }
}

\references{
  APG IV [Angiosperm Phylogeny Group] (2016) An update of the
  Angiosperm Phylogeny Group classification for the orders and
  families of flowering plants: APG IV. \emph{Bot. J. Linnean Soc.}
  \strong{181}: 1--20.

  Chase, M.W. & Reveal, J. L. (2009) A phylogenetic classification of
  the land plants to accompany APG III. \emph{Bot. J. Linnean Soc.}
  \strong{161}: 122--127.

}
\examples{
data(BCI.taxon)
}
\keyword{datasets}
