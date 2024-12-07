\name{write_broadpeaks}
\alias{write_broadpeaks}
\alias{write_seacr}

\title{
Write a BED file for MACS2 broadpeaks and SEACR peaks
}

\description{
This is a wrapper to the export functions defined in \emph{rtracklayers}
}

\usage{
write_broadpeaks(x, file)
write_seacr(x, file)
}

\arguments{
  \item{x}{A \code{GRanges} instance corresponding to either MACS2 broadPeaks or SEARC spec.}
  \item{file}{A path to a file to write \code{x}.}
}

\seealso{
See Also as \code{\link{plyranges::write_bed()}} and \code{\link{rtracklayer::BEDFILE()}}
}

\examples{

  \dontrun{

  }
}
