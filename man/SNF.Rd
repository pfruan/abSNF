\name{SNF}
\alias{SNF}
\title{
Similarity Network Fusion
}
\description{
Similarity Network Fusion takes multiple views of a network and fuses them together to construct an overall status matrix.
}
\usage{
SNF(Wall, K, t)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Wall}{
List of matrices. Each element of the list is a square, symmetric matrix that shows affinities of the data points from a certain view.
}
  \item{K}{
Number of neighbors in K-nearest neighbors part of the algorithm.
}
  \item{t}{
Number of iterations for the diffusion process.
}
}
\value{
W is the overall status matrix derived
}
\references{
B Wang, A Mezlini, F Demir, M Fiume, T Zu, M Brudno, B Haibe-Kains, A Goldenberg (2014) Similarity Network Fusion: a fast and effective method to aggregate multiple data types on a genome wide scale. Nature Methods. Online. Jan 26, 2014 \cr
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
Dist1 = dist2_w(as.matrix(data1),as.matrix(data1),weight1)^(1/2)
Dist2 = dist2_w(as.matrix(data2),as.matrix(data2),weight2)^(1/2)

# Next, construct similarity graphs
W1 = affinityMatrix(Dist1)
W2 = affinityMatrix(Dist2)

#
W = SNF(list(W1,W2), 20, 20)
}
