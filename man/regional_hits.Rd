\name{regional_hits}
\alias{regionhit_per_sample_mat}
\alias{consolidated_peak_hits}
\alias{consolidate_peaks}
\title{Read a MACS2 narrowPeak and broadPeak BED file}
\description{
  Construct a hit matrix of peaks against a collection of peaks
}

\usage{
  consolidated_peak_hits(grl, min_overlap=NULL)
  
  consolidate_peaks(grl)

  regionhit_per_sample_mat(regions, grl, min_overlap=NULL)
}

\arguments{
  \item{grl}{A list of \code{GRanges} objects or a \code{GRangesList} object representing peak ranges of multiple samples. Note that the name of \code{grl}, representing the sample ID, must be given.}
  \item{min_overlap}{The minimum amoutn of overlap between the peak and the consolidate peaks (or regions). If NULL, it will be set to the half of the size of the minimal width of the consolidated peaks (\code{grl} or \code{regions})}.
  \item{regions}{A \code{GRanges} object representing consolidated peak ranges or ranges of interest}
}

\details{
  The \code{consolidated_peak_hits} function consolidates and merges a list of peaks or a \code{GRagnesList} into a \code{GRanges} object. Suppose there are m number of samples (length of the input \code{grl}) and n number of consolidated peaks,  \code{consolidated_peak_hits} constructs a n-by-m matrix with an entry of 1 if a sample's peak hit the consolidated peaks, and 0 otherwise. What can we do with this type of matrix?
    
  The \code{consolidate_peaks} function consolidates (merges) and mergess a list of peaks or a \code{GRangesList} into a \code{GRanges} object.
  
  The \code{regionhit_per_sample_mat} function constructs a n-by-m matrix with an entry of 1 if the sample's peak hit the consolidated peaks, and 0 otherwise, where n is the length of the input regions (\code{GRanges}) and m is number of samples (length of the input list \code{grl}).
}

\value{
  A matrix with an entry of 1, if the sample’s peak hit the merged peaks, and 0 otherwise.
}


\examples{
# construct a sample_info data.frame
# 1) get SEACR bed files
library(dplyr)
library(ggplot2)
seacr_pattern = '\\.stringent.bed$'
peak_df <- data.frame(
  bed_file = list.files(system.file('extdata', package = 'peaklerrr'),
                       full.names = TRUE, pattern=seacr_pattern)) %>%
  dplyr::mutate(sample_id = stringr::str_replace(basename(bed_file), seacr_pattern, ''))      
# 2) get sample information  
sample_info <- data.frame(
  sample_id = c('chr2_Rep1_H1_CTCF', 'chr2_Rep2_H1_CTCF',
                'chr2_Rep1_H1_H3K4me3', 'chr2_Rep2_H1_H3K4me3'),
  antibody = c('CTCF', 'CTCF', 'H3K4me3', 'H3K4me3')) %>%
  dplyr::left_join(peak_df, by='sample_id')

# 3) read bed files and get a list of Granges
grl <- lapply(sample_info$bed_file, read_seacr)
names(grl) <- sample_info$sample_id

# 4) get the hit matrix against consolidated peaks of four samples
hits_mat <- consolidated_peak_hits(grl)

# 5) PCA: visualize the peaks coorelation bewteen replicates and the difference of peaks of
# bewteen antibodies, for example:
pcs <- peaklerrr:::.getPCA(hits_mat, sample_info=sample_info, n_pcs=2)
ggplot(pcs, aes(x=PC1, y=PC2, color=antibody)) +
  geom_point() +
  theme_minimal()
  

#
# another example for MACS2 narrow peaks
#
# 1) get SEACR bed files
library(dplyr)
library(ggplot2)
macs2_pattern = '\\.narrowPeak$'
sample_info <- data.frame(
  bed_file = list.files(system.file('extdata', package = 'peaklerrr'),
                       full.names = TRUE, pattern=macs2_pattern)) %>%
  dplyr::mutate(sample_id = stringr::str_replace(basename(bed_file), 
                                                 macs2_pattern, '')) %>%
  dplyr::mutate(antibody = c("CTCF", "H3K4me3", "CTCF", "H3K4me3"))                                                 
grl <- lapply(sample_info$bed_file, read_macs2_narrow)
names(grl) <- sample_info$sample_id

# 4) get the hit matrix against consolidated peaks of four samples
hits_mat <- consolidated_peak_hits(grl)

# 5) PCA: visualize the peaks coorelation bewteen replicates and the difference 
# across antibodies, for example:
pcs <- peaklerrr:::.getPCA(hits_mat, sample_info=sample_info, n_pcs=2)
ggplot(pcs, aes(x=PC1, y=PC2, color=antibody)) +
  geom_point() +
  theme_minimal()
  
}

\seealso{
  \code{\link{read_seacr}}, \code{\link{read_macs2_narrow}}, \code{\link{read_macs2_broad}}
}


\keyword{function}