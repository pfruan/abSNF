\name{post_process}
\alias{post_process}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Post process of enhanced network.
}
\description{
Denoise the enhanced network and make it symmetric.
}
\usage{
post_process(se,percent=0.9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{se}{
The enhanced disease specific network.
}
  \item{percent}{
what percentage of edges to be considered as noise.
}
}
\value{
The denoised and symmetric enhanced disease-specific network.
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
