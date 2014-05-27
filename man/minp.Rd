\name{minp}
\alias{minp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
P-values adjustment in permutation
}
\description{
Returns the p-values in each permutation.
}
\usage{
minp(p, rank, n.GeneSets, nbPerm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
input p-values.
}
  \item{rank}{
the rank of the p-values.
}
  \item{n.GeneSets}{
the number of genes in a given gene set.
}
  \item{nbPerm}{
the number of permutation times.
}
}

\value{
a permutation p-value matrix.
}
