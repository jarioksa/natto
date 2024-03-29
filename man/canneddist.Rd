% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/canneddist.R
\encoding{UTF-8}
\name{canneddist}
\alias{canneddist}
\title{Canned Dissimilarities with their Vernacular Names}
\usage{
canneddist(x, method, help = FALSE)
}
\arguments{
\item{x}{Input data.}

\item{method}{Vernacular name for a dissimilarity index.}

\item{help}{List available indices and their definitions instead of
calculating dissimilarities.}
}
\value{
Function returns a \code{\link{dist}} object of dissimilarities.
}
\description{
Function is a storehouse for dissimilarity indices that can be
called by their vernacular names. The function is based on
\code{\link[vegan]{designdist}} (\CRANpkg{vegan} package).
}
\details{
Function wraps popular dissimilarity indices for
\code{\link[vegan]{designdist}} allowing these to be called by
their popular names. It can have synonymous names for one
dissimilarity index. Use argument \code{help=TRUE} to see the
current selection of indices and their definitions.

The function uses the main notation of
\code{\link[vegan]{designdist}} where terms are based on sums and
paired minima or sum of squares and crossproducts of pairs of
sampling units. For vectors \code{x} and \code{y} the
\code{"quadratic"} terms are \code{J = sum(x*y)}, \code{A =
sum(x^2)}, \code{B = sum(y^2)}, and \code{"minimum"} terms are
\code{J = sum(pmin(x,y))}, \code{A = sum(x)} and \code{B = sum(y)},
and \code{"binary"} terms are either of these after transforming
data into binary form (number of shared species \code{J}, and
number of species for each row, \code{A, B}.). Number of columns
(species) is notated as \code{P}, and the number of sampling units
is \code{N}.

There is a huge number of indices, and the current selection is far
from comprehensive (but it can easily expanded). See References for
sources. Many sources use different notation, but they were changed
to the notation described above. For instance, in popular (but
strange) 2x2 contingency table notation for binary data \code{a =
J}, \code{b = A-J}, \code{c = B-J}, \code{d = P-A-B+J}. Some of
formulae may be surprising, but they are mathematically equivalent
to traditional ones. I challenge you to inspect Euclidean distance,
and once you see how it is derived, try Chord distance.
}
\examples{
data(spurn)
## Ochiai dissimilarity
canneddist(spurn, "ochiai")

}
\references{
\enc{Hubálek}{Hubalek}, Z. (1982). Coefficients of association and
similarity, based on binary (presence-absence) data: an
evaluation. \emph{Biological Review} 57, 669--689.

Legendre, P. & Legendre, L. (2012). \emph{Numerical Ecology.} 3rd
English Ed., Elsevier.

Yue, J.C. & Clayton, M.K. (2005). A similarity measure based on
species proportions. \emph{Communications in Statistics Theory and
Methods} 34, 2123--2131. \doi{10.1080/STA-200066418}.
}
\seealso{
Function is a wrapper to
\code{\link[vegan]{designdist}}. \CRANpkg{vegan} function
\code{\link[vegan]{betadiver}} is a similar collection of canned
indices for beta diversity, and many of these are well-known
dissimilarity indices.
}
\author{
Jari Oksanen
}
