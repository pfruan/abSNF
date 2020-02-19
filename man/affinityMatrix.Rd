\name{affinityMatrix}
\alias{affinityMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Affinity matrix calculation
}
\description{
Computes affinity matrix from a generic distance matrix
}
\usage{
affinityMatrix(Diff, K, sigma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Diff}{
Distance matrix
}
  \item{K}{
Number of nearest neighbors
}
  \item{sigma}{
Variance for local model
}
}
\value{
Returns an affinity matrix that represents the neighborhood graph of the data points.
}
\references{
B Wang, A Mezlini, F Demir, M Fiume, T Zu, M Brudno, B Haibe-Kains, A Goldenberg (2014) Similarity Network Fusion: a fast and effective method to aggregate multiple data types on a genome wide scale. Nature Methods. Online. Jan 26, 2014  \cr
Using Association Signal Annotations to boost Similarity Network Fusion (2018), Peifeng Ruan, Ya Wang, Ronglai Shen, Shuang Wang.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{

#load data
data(data1)
data(data2)
data(weight1)
data(weight2)

#standard normalization of the datasets
data1 = standardNormalization(data1)
data2 = standardNormalization(data2)

# Calculate boosted distance matrices(here we calculate Euclidean Distance, 
Dist1 = dist2_w(as.matrix(data1),as.matrix(data1),weight1)
Dist2 = dist2_w(as.matrix(data2),as.matrix(data2),weight2)

# Next, construct similarity graphs
W1 = affinityMatrix(Dist1)
W2 = affinityMatrix(Dist2)

}
