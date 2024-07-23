\name{peakable}
\alias{peakable}
\title{peakable}

\description{
  peakable provides functions facilitating CUT&RUN-seq peak caller QC and downstream analysis
}

\details{
  peakable provides three categories of functions for (1) importing 
  bed files from MACS2 and SEACR peak callers, (2) performing 
  conanical correlation analysis of peaksets, and (3) counting reads on peak ranges. Major functions includes
  \describe{
    \item{read_macs2_narrow}{Import MACS2 narrow peaks file in BED6+4 and format the intervals as \code{GRagnes}}
    \item{read_macs2_broad}{Import MACS2 broad peaks file in BED6+4 and format the intervals as \code{GRagnes}}
    \item{read_seacr}{Import SEACR's non-standard peaks BED files and format the intervals as \code{GRanges}}
    \item{find_consensus}{The \code{find_consensus_seacr} and \code{find_consensus_macs2} are wrapper functions to find consensus peaks between two SEACR or MACS2 narrow peak sets with adjusted the corresponding AUC/maximun signals and regions (SEACR) or scores/summit (MACS2).} 
    \item{consolidated_peak_hits}{Consolidates and merges a list of peak ranges or a GRagnesList into a GRanges object and constructs a n(consolidated peaks)-by-m(samples) matrix with an entry of 1 if the sample's peaks hit the consolidated peaks, and 0 otherwise.}
    \item{extract_summit}{The \code{extract_summit_macs2} and \code{extract_summit_seacr}functions extract summit information from the metadata columns of SEACR or MASC2 peaks in GRanges, respectively, and construct summit range intervals as \code{GRanges}.}
    \item{find_overlaps_venn}{A wrapper function of \code{ggVennDiagram::ggVennDiagram} and \code{plyranges::find_overlaps} to draw venn diagram of overlaps between two ragnes.}
    \item{peak_read_count}{A wrapper function of \code{GenomicAligments::summerizeOverlaps} to count reads overlap the feature ranges.}
  }
All the test peak files were profiled from GEO series GSM3391651, GSM3391652, GSM3391655, and GSM3391656.
}

\references{
  \link{https://genome.ucsc.edu/FAQ/FAQformat.html#format12}

  \link{https://github.com/FredHutch/SEACR/}
  
  \code{GenomicAlignments::\link{summerizeOverlaps}}
  
  \code{plyranges::\link{find_overlaps}}
  
  \code{ggVennDiagram::\link{ggVennDiagram}}
}

\section{Author}{
  Chao-Jen Wong <chao-jen.wong@seattlechildrens.org>
}

\keyword{package}