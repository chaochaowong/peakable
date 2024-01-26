\name{read_macs2}
\alias{read_macs2, read_macs2_narrow, read_macs2_broad}
\title{Read a MACS2 narrowPeak and broadPeak BED file}
\description{
  This is a wrapper to the `import.bed` function defined in _rtracklayer_, 
  modified to accommodate the importing of MACS2 peaks files, 
  which include extra columns beyond those found in the standard BED format. 
  
  Read intervals and format as `GRanges`.
}
\usage{
  read_macs2_narrow(file, ...)
  
  read_macs2_broad(file, ...)
}
\arguments{
  \item{file}{A path to a file or a connection.}
  \item{...}{Extra arguments pass to rtracklayer::import.bed().}
  ...
}
\details{
  This is a wrapper to the `import.bed` function defined in _rtracklayer_, 
  modified to accommodate the importing of MACS2 peaks files, 
  which include extra columns beyond those found in the standard BED format. 
  See https://genome.ucsc.edu/FAQ/FAQformat.html#format12 for details.
  
  The `read_macs2_narrow` function parses the MACS2 narrowPeak BED format.
  
  The `read_macs2_board` function parses the MACS2 broadPeak BED format.
}
\value{
  A GRanges object
}


\examples{
# narrowPeaks
narrow_file <- system.file('extdata', 'chr2_Rep1_H1_CTCF_peaks.narrowPeak', 
                            package='peaklerrr')
gr <- read_macs2_narrow(narrow_file)
gr

# broadPeaks
broad_file <- system.file('extdata, 'chr2_Rep1_H1_H3K27me3_peaks.broadPeak',
                           package='peaklerrr')
gr <- read_macs2_broad(broad_file)
gr

}
\references{
  https://genome.ucsc.edu/FAQ/FAQformat.html#format12 
}

\seealso{
  \code{\link{rtracklayer::import.bed}}
}

\keyword{function}

