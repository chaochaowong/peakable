\name{read_seacr}
\alias{read_seacr}
\title{Read a SEACR BED file}
\description{
  Import SEACR's non-standard peaks BED files and format intervals as \code{GRanges}
}

\usage{
  read_seacr(file)
}

\arguments{
  \item{file}{A path to a file or a connection.}
}

\details{
  The resulting peak files by SEACR is not in the standard BED6 format. \code{read_seacr} parses the SEACR peak format and formats the ranges into GRanges along with columns, including "AUC", "max.signal", and "max.signal.region", in addition to "num" if provided. See https://github.com/FredHutch/SEACR/ for details. 
}

\value{
  A GRanges object
}


\examples{
# SEACR peak BED file
seacr_file <- system.file('extdata', 
                          'chr2_Rep1_H1_CTCF.stringent.bed', 
                            package='peaklerrr')
gr <- read_seacr(seacr_file)
gr

}

\references{
  https://github.com/FredHutch/SEACR/
}

\keyword{function}

