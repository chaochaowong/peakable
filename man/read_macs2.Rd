\name{read_macs2}
\alias{read_macs2}
\alias{read_macs2_broad}
\alias{read_macs2_narrow}
\title{Read a MACS2 narrowPeak and broadPeak BED file}
\description{
  This is a wrapper to the \code{import.bed} function defined in \pkg{rtracklayer}, modified to accommodate the importing of MACS2 peaks file in BED6+4. Intervals are formatted as \code{GRanges}.
}

\usage{
  read_macs2_narrow(file,
                    drop_chrM = FALSE,
                    keep_standard_chrom = FALSE,
                    species = NULL, ...)

  read_macs2_broad(file,
                   drop_chrM = FALSE,
                   keep_standard_chrom = FALSE,
                   species = NULL, ...)
}

\arguments{
  \item{file}{A path to a file or a connection.}
    \item{drop_chrM}{A logical value indicating whether to drop peaks of chromosome M. Default: FALSE.}
  \item{keep_standard_chrom}{A logical value indicating whether to keep only standard seqlevels. Default: FALSE}
  \item{species}{The species of the organism. Supported sepcies can bee seen with \code{names(genomeStyles())}.}
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
                            package='peakable')
gr <- read_macs2_narrow(narrow_file)
gr

# test for the extra setting; won't change the peak size
gr <- read_macs2_narrow(narrow_file,
                        drop_chrM = TRUE,
                        keep_standard_chrom = TRUE,
                        species = 'Homo_sapiens')
gr

# MACS2 broadPeaks
broad_file <- system.file('extdata',
                          'chr2_Rep1_H1_H3K27me3_peaks.broadPeak',
                           package='peakable')
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
