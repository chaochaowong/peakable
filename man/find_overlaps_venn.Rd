\name{find_overlaps_venn}
\alias{find_overlaps_venn}
\title{find_overlaps_venn}
\description{
  A wrapper function of \code{ggVennDiagram::ggVennDiagram} and \code{plyranges::find_overlaps} to draw venn diagram of overlaps between two ranges
}

\usage{
  find_overlaps_venn(x, y, 
                     label_x = NULL, label_y = NULL,
                     maxgap = -1L, minoverlap = 1L, 
                     ...)
}

\arguments{
  \item{x, y}{\code{GRanges} objects representing ranges}
  \item{label_x, label_y}{Labels for two ranges. Default to x and y.}
  \item{maxgap}{The maximum gap between two ranges}
  \item{minoverlap}{The minimum amout of overlap between two ranges}
  \item{...}{Arguements pass to \code{ggVennDiagram::\link{ggVennDiagram}}}
}

\value{
  A \code{ggplot} object
}

\examples{
# SEACR peaks
library(ggplot2)
seacr_file_x <- 
  system.file('extdata',
              'chr2_Rep1_H1_CTCF.stringent.bed',
              package='peakable')
seacr_file_y <- 
  system.file('extdata',
              'chr2_Rep2_H1_CTCF.stringent.bed',
              package='peakable')
                            
x <- peakable::read_seacr(seacr_file_x)
y <- peakable::read_seacr(seacr_file_y)

minoverlap <- min(min(width(x)), min(width(y))) / 2
find_overlaps_venn(x, y, minoverlap = minoverlap)

}

\references{
  \link{https://genome.ucsc.edu/FAQ/FAQformat.html#format12 }
  
  \link{https://github.com/FredHutch/SEACR/}
}




\keyword{function}