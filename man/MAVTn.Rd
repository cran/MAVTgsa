\name{MAVTn}
\alias{MAVTn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
OLS, Hottelling's T2 and MANOVA with n contrasts
}
\description{
             A gene set analysis functions for computiong the p-values for one-sided test (OLS) and two-sided test (multivariate analysis of variance).
             If the experimental conditions are equal to 2, the p-value for Hotelling's t square test is calculated.
             If the experimental conditions are great than 2, the p-value for Wilks' Lambda is deterimined and post-hoc test is reported too.
             The p-value for individual gene test of significant gene sets are also listed.
}
\usage{
MAVTn(DATA, GS, MCP = 1, alpha = 0.01, nbPerm = 5000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DATA}{
an (m+1) x n gene expression data matrix with n samples in columns. The first row contains the information of experimental condition of each sample.  The genes are ex-pressed in the rest m rows.
}
  \item{GS}{
an m x k binary matrix with code (0, 1), where k is the number of gene sets.  Each column represents a pre-defined gene set.
}
  \item{MCP}{
the choice for one of three multiple comparison methods, Dunnett = 1, Tuckey = 2, Sequential pairwise = 3.
}
  \item{alpha}{
the significant level
}
  \item{nbPerm}{
the number of permutation specified
}
}

\value{
The p-values of OLS and MANOVA test are returned. If there is any significant gene set, the p values for individual genes in the gene set will be reported.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{

#------------simulate data matrix---------------#

data(data)
data(GS)

MAVTn(data,GS,MCP=1, nbPerm = 100)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MANOVA test}

