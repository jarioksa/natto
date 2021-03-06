% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diverclust.R
\name{diverclust}
\alias{diverclust}
\title{beta and gamma Diversity Clustering}
\usage{
diverclust(
  x,
  renyi = 1,
  equalize = TRUE,
  beta = TRUE,
  hill = FALSE,
  trace = TRUE
)
}
\arguments{
\item{x}{Community data.}

\item{renyi}{The scale of Renyi diversity index as defined in
\code{\link[vegan]{renyi}}. Value \code{0} gives diversity based on
species richness, \code{1} diversity based on Shannon entropy,
\code{2} diversity based on Simpson index, and \code{Inf} diversity
based on Berger-Parker index. Other non-negative values are also
allowed, but they may not define a well known standard diversity
index.}

\item{equalize}{Equalize data rows so that they can be meaningfully
pooled and averaged. If this is \code{FALSE}, raw data will be
used, and beta diversities may be negative. The equalization
depends on the value of \code{renyi}, and each row of \code{x} is
divided by \code{(rowSums(x^renyi))^(1/renyi)}. In borderline cases
\code{renyi=0} the data are presence-absence tranformed (but can be
left untransformed), and with \code{renyi=Inf} the rows are divided
by row maxima (see \code{\link[vegan]{decostand}}).}

\item{beta}{Use beta diversities: the average alpha diversity of
cluster members is subtracted from the pooled diversity. If this is
\code{FALSE},the clustering is based on pooled diversity, also
known as gamma diversity.}

\item{hill}{Use Hill numbers instead of Renyi diversity (see
\code{\link[vegan]{renyi}}). For \code{renyi = 0} these are species
richess values instead of their logarithms, and for \code{renyi =
1} they are exponents of Shannon diversity instead of Shannon
diversities. In general, a Hill number is an exponent of Renyi
diversity. The Hill numbers may not be strictly additive, and beta
diversities may be negative (except with \code{renyi = 0}).}

\item{trace}{Trace calculations. Either logical or an integer:
\code{trace=2} also traces merges.}
}
\value{
Function returns an \code{\link{hclust}} object.
}
\description{
Function forms hierarchic clustering so that each fusion minimizes
beta deiversity or gamma diversity.
}
\details{
The function forms clusters so that pooled diversity or
its change to the baseline alpha diversity is minimized at each
level of clustering. The change of pooled diversity with respect to
the baseline alpha diversity is called beta diversity, and the
overall diversity is called gamma diversity. For beta diversity,
the clustering implies an additive partitioning.
}
\seealso{
\code{\link{hclust}} for cluster analysis and its support
methods, \code{\link[vegan]{renyi}} for estimating Renyi
diversities and Hill numbers, and \code{\link[vegan]{adipart}} for
related additive partitioning of beta diversity.
}
\author{
Jari Oksanen.
}
