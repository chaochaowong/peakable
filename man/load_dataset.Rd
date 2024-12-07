\name{load_dataset}
\alias{load_dataset}

\title{
Load a sample dataset that is used for the vignette and examples in the man pages
}
\description{
Construct a list instance containing a \code{list} of peak ranges in \code{GRanges} and a \code{data.frame} representing the metadata of the peakset.
}

\usage{
load_dataset()
}


\value{
A list of the following items:
  \item{df}{A \code{data.frame} of sample and peak-calling (MACS2 narrow) information including the sample id, cell line, anitoby, peak number, peak caller, and path to peak bed files.}
  \item{grl}{A list of peak ranges (\code{GRanges}) corresponding to the samples in \code{df}.}

}

\author{
Chao-Jen Wong
}

\examples{
peaksets <- load_dataset()
peaksets$df
peaksets$grl
}
