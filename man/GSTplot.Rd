\name{GSTplot}
\alias{GSTplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GST plot
}
\description{
The GST plot displays the relative direction (in two conditions) and statistics ranking for genes in a gene set.
}
\usage{
GSTplot(data, gs, geneset.name = NULL, alpha = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a gene expression data matrix with samples in columns
}
  \item{gs}{
a binary matrix coded 0 or 1 with genes in rows
}
  \item{geneset.name}{
The name of the given gene set
}
  \item{alpha}{
the significance level
}
}


\references{
Chen,J.J. et al. (2007) Significance analysis of group of genes in expression profiling studies, Bioinformatics, 23, 2104.

Tsai,C.A. et al. (2009) Multivariate analysis of variance test for gene set analysis. Bioinformatics, 25, 897. 
}
\author{
Chih-Yi Chien, Chen-An Tsai, Ching-Wei Chang, and James J. Chen
}
\note{
R > 2.13.2
}

