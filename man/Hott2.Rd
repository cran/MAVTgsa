\name{Hott2}
\alias{Hott2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hottelling's T square
}
\description{
To compute Hotelling's T square statistic for multivariate analysis of variance using Shrinkage covariance matrix estimates.
}
\usage{
Hott2(x, y, var.equal = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Data matrix; row is sample; each column is variable(gene) 
}
  \item{y}{
Vector defining two-group of the samples.
}
  \item{var.equal}{
Logical.
}
}

\value{
Hotelling's T square statistic is calculated.
}
\references{
Chen,J.J. et al. (2007) Significance analysis of group of genes in expression profiling studies, Bioinformatics, 23, 2104.

Tsai,C.A. et al. (2009) Multivariate analysis of variance test for gene set analysis. Bioinformatics, 25, 897. 
}
\author{
Chen-An Tsai, James J. Chen, Ching-Wei Chang, and Chih-Yi Chien
}
\note{
R > 2.13.2
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.


