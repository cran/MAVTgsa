\name{MAVTp}
\alias{MAVTp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Random Forests-based procedure
}
\description{
 A Random Forests-based procedure is to identify gene sets that can accurately predict samples from different experimental conditions or are associated with the continuous phenotypes. 
}
\usage{
MAVTp(DATA, GS, nbPerm = 5000, numoftree = 500, type = c("cont", "cate"), impt = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DATA}{
a gene expression data matrix with samples in columns. The first row contains the information of the experimental condition of each sample.  The remaining rows contain gene expression. 
}
  \item{GS}{
an m x k binary matrix with code (0, 1), where k is the number of gene sets.  Each column represents a pre-defined gene set.
}
  \item{nbPerm}{
the number of permutation specified
}
  \item{numoftree}{
the number of trees to grow
}
  \item{type}{
 This can be one of "cont" (continuous phenotypes) and "cate" (categorical phenotypes).
}
  \item{impt}{
 If TRUE (default), the importance measurement will be output. 
}
}

\value{
A list of the p-values of random forests for GSA. The importance measurement of individual genes for those significant gene sets will also be output when impt is set TRUE.
}
\references{
H.M. Hsueh, et al. (2013) Random forests-based differential analysis of gene sets for gene expression data. Gene, 518, 179-186.}

\author{
Chih-Yi Chien, Chen-An Tsai, Ching-Wei Chang, and James J. Chen 
}
\note{
R > 2.14.0
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{

data(data)
data(GS)
a=proc.time()
MAVTp(data,GS, nbPerm = 50, numoftree = 500, type = "cate", impt = TRUE)
proc.time()-a

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Random forests}

