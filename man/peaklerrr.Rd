\name{peaklerrr}
\alias{peaklerrr}
\title{peaklerrr}

\description{
  peaklerrr provides functions facilitating CUT&RUN-seq peak caller QC 
}

\details{
  peaklerrr provides two categories of functions for (1) importing 
  bed files from MACS2 and SEACR callers, and (2) performing 
  conanical correlation analysis of peaksets. Functions includes
  \itemize{
    \item{read_macs2_narrow`{Import MACS2 narrow peaks file in BED6+4 and format the intervals as \code{GRagnes}}
    \item{read_macs2_broad}{Import MACS2 broad peaks file in BED6+4 and format the intervals as \code{GRagnes}}
    \item{read_seacr}{Import SEACR's non-standard peaks BED files and format the intervals as \code{GRanges}}
    \item{find_consensus_seacr}{Find consensus peaks between two peak sets with adjust score, summit (MACS2), or maximun signals and regions (SEACR)} 
    \item{consolidated_peak_hits}{Consolidates and merges a list of peaks or a GRagnesList into a GRanges object and constructs a n(consolidated peaks)-by-m(samples) matrix with an entry of 1 if the sample's peaks hit the consolidated peaks, and 0 otherwise.}
    \item{extract_summit}{The \code{extract_summit_macs2} and \code{extract_summit_seacr} functions extract summit information from the metadata columns of SEACR or MASC2 peaks in GRanges, respectively, and construct summit range intervals as \code{GRanges}.}
    \item{find_overlaps_venn}{}

  }
  
  All the test peak files were profiled from GEO series GSM3391651, GSM3391651, 
  and GSM3391656.
}

\section{Author}{
  \itemize{
    \item{Chao-Jen Wong (chaojen.wong@seattlechildrens.org)}
  }
}

\keyword{package}