\name{cos_similarity_by}
\alias{cos_similarity_by}
\title{Cosine similarity between defined replicates}
\description{
  A wraper function of \link{\code{cos_similarity}}. Take the list of GRanges, data.frame, and variables to group by, \code{cos_similarity_by} performs \code{cos_similarity} "by group".
}
\usage{
cos_similarity_by(sample_df,
                  peaks_grl,
                  sim_group_by)
}
\arguments{
  \item{sample_df}{A data.frame or extension (e.g. \code{tibble}) representing the metadata of the list of peak ranges }
  \item{peaks_grl}{A list of \code{GRanges} of peak ranges}
  \item{sim_group_by}{variables (e.g. column names in \code{sample_df}) to group the peak ranges}
}
\value{
  The cosine similarity between peaks ranges that are "group by" \code{sim_group_by}.
}
\examples{
peaksets <- load_dataset()
# cosine similarity between replicates defined by "cell_line" and
# "antibody" colums
cos_sim <- cos_similarity_by(sample_df = peaksets$df,
                             peaks_grl = peaksets$grl,
                             sim_group_by = c('cell_line', 'antibody'))
cos_sim
}

