\name{peak_read_count}
\alias{peak_read_count}
\title{Perform counting paired-end mapped read on peak ranges}
\description{
  A wrapper function of \code{GenomicAlignments::summarizeOverlaps} for faciliting read counting for peak ranges. 
}

\usage{
  peak_read_count(features, 
                  sample_df, 
                  spike_in_norm = FALSE)
}

\arguments{
  \item{features}{A \code{GRanges} object of peak regions of interest.}
  \item{sample_df}{A \code{data.frame} object obtaining sample information: column sample_id", and "bam_file" are mandatary.}
  \item{spike_in_norm}{A logical value to reinforce spike-in normaliztion on read counts. If TRUE, sample_df must obtain a "spike_in_factor" column.}
}

\details{
  The \code{peak_read_count} is a wrapper function of \code{summarizeOverlaps} for facilitating read counting for peak ranges of interest from pair-ended bam files. The sample_df must obtain \code{sample_id} and \code{bam_file}. Note that the \code{summarizeOverlaps} parameters were set to take pair-ended reads only and ignore strand.
}

\value{
  A \code{RangedSummarizedExperiment} object with the "sample_info" appended to the column data, \code{colData}.
}


\examples{
  \donotrun{
# sample info
bam_pattern <- '\\.bam$'
bam_df <- data.frame(
  bed_file = list.files(system.file('extdata', package = 'peaklerrr'),
                       full.names = TRUE, pattern='\\.bam$')) %>%
  dplyr::mutate(sample_id = 
                stringr::str_replace(basename(bam_file), 
                bam_pattern, ''))
                
seacr_pattern = '\\.stringent.bed$'                
bed_df <- peak_df <- data.frame(
  bed_file = list.files(system.file('extdata', package = 'peaklerrr'),
                       full.names = TRUE, pattern=seacr_pattern)) %>%
  dplyr::mutate(sample_id = 
                stringr::str_replace(basename(bed_file), 
                seacr_pattern, ''))         
                
# 2) get sample information

sample_info <- data.frame(
  sample_id = c('chr2_Rep1_H1_CTCF', 'chr2_Rep2_H1_CTCF',
                'chr2_Rep1_H1_H3K4me3', 'chr2_Rep2_H1_H3K4me3'),
  antibody = c('CTCF', 'CTCF', 'H3K4me3', 'H3K4me3')) %>%
  dplyr::left_join(bam_df, by='sample_id') %>%
  dplyr::left_join(bed_df, by='sample_id')
  
# define feature: 

x <- peaklerrr::read_seacr(sample_info$bed_file[1])
y <- peaklerrr::read_seacr(sample_info$bed_file[3]) 
features <- GenomicRanges::reduce(x, y)

# count reads for peak ranges              
se <- peak_read_count(features, sample_df, 
                      spike_in_norm = FALSE)
se            

# visualize the correlation: Pearson and PCA

  }
}

\references{
  \link{https://github.com/FredHutch/SEACR/}
}

\keyword{function}

