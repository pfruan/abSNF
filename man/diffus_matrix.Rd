\name{diffus_matrix}
\alias{diffus_matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Network enhancement by diffusion process.
}
\description{
Enhances a comprehensive gene network for a disease using a gene-gene similarity network derived from a type of omics data of the disease that reflects correlations/interactions between genes for the disease.
}
\usage{
diffus_matrix(s0, adjacency, alpha=0.75,iter=10, difference=1e-6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s0}{
The original general gene network.
}
  \item{adjacency}{
A gene-gene similarity network derived from a type of omics data of a disease.
}
  \item{alpha}{
A regularization parameter representing weights for signal sources where α = 0 means no disease-specific enhancement.
}
  \item{iter}{
Number of iterations.
}
  \item{difference}{
A parameter for converenge.
}
}
\value{
The enhanced disease-specific network.
}
\references{
DiSNEP: a Disease-Specific gene Network Enhancement to improve Prioritizing candidate disease genes (2020), Peifeng Ruan, Shuang Wang.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{

#enhance the general network s0 into a disease specific network by diffusion on similarity network generated from disease omics data.
se=diffus_matrix(s0,adjacency)

#denoise the enhanced network and make it symmetric.
se_post=post_process(se)

#prioritize the disease association signals.
p1=diffus_vec(signals,se_post,"pvalue")

}
