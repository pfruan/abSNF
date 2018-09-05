\name{dist2}
\alias{dist2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pairwise Euclidean distances
}
\description{
Computes the Euclidean distances between all pairs of data point given
}
\usage{
dist2(X, C)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A data matrix where each row is a different data point
}
  \item{C}{
A data matrix where each row is a different data point. If this matrix is the same as X, pairwise distances for all data points are computed.
}
}
\value{
Returns an N x M matrix where N is the number of rows in X and M is the number of rows in M. element (n,m) is the squared Euclidean distance between nth data point in X and mth data point in C
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

#load data
data(data1)
data(data2)

# Calculate distance matrices(here we calculate Euclidean Distance, 
Dist1 = dist2(as.matrix(data1), as.matrix(data1))
Dist2 = dist2(as.matrix(data2), as.matrix(data2))

}
