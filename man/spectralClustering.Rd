\name{spectralClustering}
\alias{spectralClustering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spectral Clustering
}
\description{
Perform the famous spectral clustering algorithms. There are three variants. The default one is the third type. 
}
\usage{
spectralClustering(affinity, K, type = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{affinity}{
Similarity matrix
}
  \item{K}{
Number of clusters
}
  \item{type}{
The variants of spectral clustering to use.
}
}
\value{
A vector consisting of cluster labels of each sample.
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
Dist1 = dist2_w(as.matrix(data1),as.matrix(data1),weight1)^(1/2)
Dist2 = dist2_w(as.matrix(data2),as.matrix(data2),weight2)^(1/2)

# Next, construct similarity graphs
W1 = affinityMatrix(Dist1)
W2 = affinityMatrix(Dist2)

#
W = SNF(list(W1,W2), 20, 20)

# 
labels = spectralClustering(W, 3)
}
