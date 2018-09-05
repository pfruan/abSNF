\name{dist2_w}
\alias{dist2_w}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Boosted Pairwise Euclidean distances
}
\description{
Computes the boosted Euclidean distances between all pairs of data point given
}
\usage{
dist2_w(X, C, weight)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A data matrix where each row is a different data point
}
  \item{C}{
A data matrix where each row is a different data point. If this matrix is the same as X, pairwise distances for all data points are computed.
}
  \item{weight}{
A sequence indicating the weights of each features.
}
}
\value{
Returns an N x M matrix where N is the number of rows in X and M is the number of rows in M. element (n,m) is the squared Euclidean distance between nth data point in X and mth data point in C
}
\references{
Weighed network fusion through association signal annotations, Peifeng Ruan, Ya Wang, Ronglai Shen, Shuang Wang.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

#load data
data(data1)
data(data2)
data(weight1)
data(weight2)

# Calculate boosted distance matrices(here we calculate Euclidean Distance, 
Dist1 = dist2_w(as.matrix(data1), as.matrix(data1),weight1)
Dist2 = dist2_w(as.matrix(data2), as.matrix(data2),weight2)

}
