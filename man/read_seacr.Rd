\name{read_seacr}
\alias{read_seacr}
\title{Read a SEACR BED file}
\description{
  Import SEACR's non-standard peaks BED files and format intervals as \code{GRanges}
}

\usage{
  read_seacr(file, 
             drop_chrM = FALSE, 
             keep_standard_chrom = FALSE,
             species = NULL)
}

\arguments{
  \item{file}{A path to a file or a connection.}
  \item{drop_chrM}{A logical value indicating whether to drop peaks of chromosome M. Default: FALSE.}
  \item{keep_standard_chrom}{A logical value indicating whether to keep only standard seqlevels. Default: FALSE}
  \item{species}{The species of the organism. Supported sepcies can bee seen with \code{names(genomeStyles())}.}
}

\details{
  The resulting peak files by SEACR is not in the standard BED6 format. \code{read_seacr} parses the SEACR peak format and formats the ranges into GRanges along with columns, including "AUC", "max.signal", and "max.signal.region", in addition to "num" if provided. See \link{https://github.com/FredHutch/SEACR/} for details. 
}

\value{
  A GRanges object
}


\examples{
# SEACR peak BED file
seacr_file <- system.file('extdata', 
                          'chr2_Rep1_H1_CTCF.stringent.bed', 
                          package='peakable')
gr <- read_seacr(seacr_file)
gr

# drop seqlevel chrM and only keep standard seqlevels
gr <- read_seacr(seacr_file, 
                 drop_chrM = TRUE, 
                 keep_standard_chrom = TRUE,
                 species = 'Homo_sapiens')
gr                 

}

\references{
  \link{https://github.com/FredHutch/SEACR/}
}

\keyword{function}

