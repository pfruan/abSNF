\name{s0}
\alias{s0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A example of a general gene network.
}
\description{
A example of a general gene network with 1000 genes randomly selected from STRING network.
}


\references{
DiSNEP: a Disease-Specific gene Network Enhancement to improve Prioritizing candidate disease genes (2020), Peifeng Ruan, Shuang Wang.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{

#enhance the general network s0 into a disease specific network by diffusion on similarity network generated from disease omics data.
se=diffus_matrix(s0,adjacency)

}
