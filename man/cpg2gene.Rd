\name{cpg2gene}
\alias{cpg2gene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mapping CpG signals to genes.
}
\description{
Mapping association signals on CpG sites to the correponding genes.
}
\usage{
cpg2gene(cpg, method='smallest')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cpg}{
A data frame for the CpG sites' signals with three columns: "gene", "cpg" and "pvalue"
}
  \item{method}{
Use which method to map the ssociation signals on CpG sites to the correponding genes.
}
}
\value{
The data frame with one column of gene name and one column of their p values.
}
\references{
DiSNEP: a Disease-Specific gene Network Enhancement to improve Prioritizing candidate disease genes (2020), Peifeng Ruan, Shuang Wang.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{

#gene_signals=cpg2gene(cpg)

}
