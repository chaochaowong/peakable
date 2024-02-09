\name{find_consensus}
\alias{find_consensus_seacr}
\alias{find_consensus_macs2}

\title{Find consensus ranges between two peak ranges}
\description{
  A wrapper of \code{plyranges::find_overlaps} to find overlaps peaks between two SEACR (or MACS2 narrow) peak ranges and carry on the modified metadata columns 
}

\usage{
  find_consensus_seacr(x, y, minoverlap = NULL) 
  
  find_consensus_macs2(x, y, minoverlap = NULL)
}

\arguments{
  \item{x, y}{\code{GRanges} objects representing two peak ranges}
  \item{minoverlap}{Default to half of the minimal width of x and y}
}

\details{
  
  The \code{find_consensus_seacr} find overlaps between to SEACR peak ranges and carry on the metadata columns corresponding to the ranges with hightest signal, (\code{max.signal}), between x and y. The metadata columns include "AUC", "max.signal"" and "max.sig.regions".
  
  The \code{find_consensus_macs2} find oerlaps between to MACS2 narrow peak ranges and carry on the metadata columns corresponding to ranges of x. (Need improvement here).
}

\value{
  A \code{GRanges} object with metadata columns (\code{mcol}) 
}


\examples{
# two SEACR peak ranges
seacr_file_x <- 
  system.file('extdata',
              'chr2_Rep1_H1_CTCF.stringent.bed',
              package='peaklerrr')
seacr_file_y <- 
  system.file('extdata',
              'chr2_Rep2_H1_CTCF.stringent.bed',
              package='peaklerrr')
                            
x <- peaklerrr::read_seacr(seacr_file_x)
y <- peaklerrr::read_seacr(seacr_file_y)
consensus <- find_consensus_seacr(x, y, minoverlap = 40L)
consensus
 
consensus <- find_consensus_seacr(x, y)
consensus

# MACS2 narrow peak ranges
macs2_file_x <- 
  system.file('extdata',
              'chr2_Rep1_H1_CTCF_peaks.narrowPeak',
               package='peaklerrr')
macs2_file_y <- 
  system.file('extdata',
              'chr2_Rep2_H1_CTCF_peaks.narrowPeak',
                package='peaklerrr')
                            
x <- peaklerrr::read_macs2_narrow(seacr_file_x)
y <- peaklerrr::read_macs2_narrow(seacr_file_y)
consensus <- find_consensus_macs2(x, y, minoverlap = 40L)
consensus

# use default minoverlap = half of the minimal width of x and y 
consensus <- find_consensus_macs2(x, y)
consensus

}

\references{
  \link{https://genome.ucsc.edu/FAQ/FAQformat.html#format12 }
  
  \link{https://github.com/FredHutch/SEACR/}
}




\keyword{function}