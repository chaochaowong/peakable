\name{extract_summmit}
\alias{extract_summit}
\alias{extract_summit_macs2}
\alias{extract_summit_seacr}
\title{Exact MACS2 narrow or SEACR peaks's summits and convert to a GRanges object}
\description{
  The \code{extract_summit_macs2} and \code{extract_summit_seacr} functions extract summit information from the metadata columns of SEACR or MASC2 peaks in \code{GRanges}, respectively, and construct summit range intervals as a \code{GRanges} object.
}

\usage{
  extract_summit_macs2(gr, summit_wid=NULL)
  
  extract_summit_seacr(gr, summit_wid=NULL)
}

\arguments{
  \item{gr}{A \code{GRanges} object representing MACS2 or SEACR peaks.}
  \item{summit_wid}{NULL (default) or a positive integer indicating the width of summit. If \code{SEACR} peaks, default to the width of 'max.signal.region'. If \code{MACS2} narrow peaks, default to 1L.}
}

\details{
  
  The \code{extract_summit_macs2} function extracts the "peak" column from the MACS2 narrowPeak \code{GRanges} and creates summit range interval as \code{GRanges}, while preserving all other columns. The width is defaulted to 1L.
  
  The \code{extract_summit_seacr} function extracts the "max.signal.region" from the SEACR peak \code{GRanges} and creates summit range intervals as \code{GRanges}, with metadata colmns "name", for arbitrary peak names, and the original "AUC" and "max.signal". The width of the "max.signal.region" is preserved, unless \code{summit_wid} is specified.
  }

\value{
  A GRanges object
}


\examples{
# MACS2 narrowPeaks
narrow_file <- 
  system.file('extdata', 
              'chr2_Rep1_H1_CTCF_peaks.narrowPeak', 
               package='peaklerrr')
gr <- read_macs2_narrow(narrow_file)
summit_macs2 <- extract_summit_macs2(gr)
summit_macs2

# SEACR peaks
seacr_file <- 
  system.file('extdata', 
              'chr2_Rep1_H1_CTCF.stringent.bed', 
              package='peaklerrr')
gr <- read_seacr(seacr_file)
summit_seacr <- extract_summit_seacr(gr)
summit_seacr
}

\references{
  \link{https://genome.ucsc.edu/FAQ/FAQformat.html#format12 }
  
  \link{https://github.com/FredHutch/SEACR/}
}




\keyword{function}