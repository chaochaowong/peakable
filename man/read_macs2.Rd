\name{read_macs2}
\alias{read_macs2}
\alias{read_macs2_broad}
\alias{read_macs2_narrow}
\title{Read a MACS2 narrowPeak and broadPeak BED file}
\description{
  This is a wrapper to the \code{import.bed} function defined in \pkg{rtracklayer}, modified to accommodate the importing of MACS2 peaks file in BED6+4. Intervals are formatted as \code{GRanges}.
}

\usage{
  read_macs2_narrow(file, ...)
  
  read_macs2_broad(file, ...)
}

\arguments{
  \item{file}{A path to a file or a connection.}
  \item{...}{Extra arguments pass to \code{rtracklayer::\link{import.bed}}.}
}

\details{
  This is a wrapper to the \code{import.bed} function defined in \pkg{rtracklayer}, modified to accommodate the importing of MACS2 peaks files in BED6+4 format. See https://genome.ucsc.edu/FAQ/FAQformat.html#format12 for details.
  
  The \code{read_macs2_narrow} function parses the MACS2 narrowPeak BED format.
  
  The \code{read_macs2_board} function parses the MACS2 broadPeak BED format.
}

\value{
  A GRanges object
}


\examples{
# MACS2 narrowPeaks
narrow_file <- system.file('extdata', 
                           'chr2_Rep1_H1_CTCF_peaks.narrowPeak', 
                            package='peaklerrr')
gr <- read_macs2_narrow(narrow_file)
gr

# MACS2 broadPeaks
broad_file <- system.file('extdata', 
                          'chr2_Rep1_H1_H3K27me3_peaks.broadPeak',
                           package='peaklerrr')
gr <- read_macs2_broad(broad_file)
gr

}

\references{
  https://genome.ucsc.edu/FAQ/FAQformat.html#format12 
}

\seealso{
  \code{rtacklayer::\link{import.bed}}
}


\keyword{function}