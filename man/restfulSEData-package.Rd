\name{restfulSEData-package}

\alias{restfulSEData-package}
\alias{restfulSEData}

\title{Example metadata for the "restfulSE" R package}

\description{
  Metadata RangedSummarizedExperiment shell for use with "restfulSE" R package 
  is available in \code{ExperimentHub}
}


\examples{
  library(ExperimentHub)
  ehub <- ExperimentHub()
  myfiles <- query(ehub, "restfulSEData")
  myfiles[[1]] #load the first resource in the list
  myfiles[["EH551"]] #load by EH id
}
