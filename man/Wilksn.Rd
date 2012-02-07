\name{Wilksn}
\alias{Wilksn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wilk's Lambda for n-group multiple comparisons
}
\description{
To compute Wilk's Lambda statistic for multivariate analysis of variance and multiple comparisons. 
}
\usage{
Wilksn(Y, class, type = c("Tukey", "Dunnett", "Sequence"), base = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
Data matrix; row is sample; each column is variable(gene). 
}
  \item{class}{
Vector defining the clinical outcome of the samples.
}
  \item{type}{
Type of contrast
}
  \item{base}{
An integer to denote which group is considered the baseline group for Dunnett contrasts.
}
}

\value{
Wilk's Lambdas for MANOVA and multiple comparisons are returned.
}
\references{
Frank Bretz, Torsten Hothorn and Peter Westfall (2010), Multiple Comparison Using R, CRC Press, Boca Raton 

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
\keyword{ Wilks' Lambda }

