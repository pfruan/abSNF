\name{standardNormalization}
\alias{standardNormalization}
\title{
Standard Normalization
}
\description{
Normalize each column of the input data to have mean 0 and standard deviation 1.
}
\usage{
standardNormalization(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The unnormalized data.
}
}
\value{
The data normalized. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

#
data(data1)

data1 = standardNormalization(data1)
}
