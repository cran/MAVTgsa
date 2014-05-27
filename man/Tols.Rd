\name{Tols}
\alias{Tols}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ordinary Least Square test
}
\description{
To compute OLS statistic for one-sided test
}
\usage{
Tols(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Data matrix; row is sample; each column is variable(gene).
}
  \item{y}{
Vector defining the clinical outcome of the samples.
}
}

\value{
Returns OLS test statistic for gene set analysis
}
\references{
Chen,J.J. et al. (2007) Significance analysis of group of genes in expression profiling studies, Bioinformatics, 23, 2104.
}
\author{
Chih-Yi Chien, Chen-An Tsai, Ching-Wei Chang, and James J. Chen 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ OLS }

