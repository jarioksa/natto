\name{BCI.taxon}
\alias{BCI.taxon}
\docType{data}
\title{
   Botanical Classification of Barro Colorado Island Species
}
\description{
  Classification of \code{\link[vegan]{BCI}} species based on
  APG III (2009) and Soltis \emph{et al.} (2011).
}
\usage{data("BCI.taxon")}
\format{
  A data frame with 225 observations on the following 7 variables.
  \describe{
    \item{\code{genus}}{Genus}
    \item{\code{family}}{APG III family}
    \item{\code{order}}{APG III order}
    \item{\code{unrank1}}{unranked level above previous}
    \item{\code{unrank2}}{unranked level above previous}
    \item{\code{unrank3}}{unranked level above previous}
    \item{\code{subclass}}{Either \code{eudicotyledons},
      \code{Magnoliidae} or \code{Monocotyledons} which sometimes are
      regarded as subclasses (but this may be too high a status). }
  }
}

\references{
  APG III [Angiosperm Phylogeny Group] (2009) An update of the Angiosperm
  Phylogeny Group classification for the orders and families of flowering
  plants: APG III. \emph{Bot. J. Linnean Soc.} \strong{161}: 105--121.

  Soltis, D.E \emph{et al.} (2011) Angiosperm phylogeny: 17 genes, 640
  taxa. \emph{Amer. J. Bot.} \strong{98}: 704--730.
}
\examples{
data(BCI.taxon)
}
\keyword{datasets}
