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
    \item{`read_macs2_narrow`}
    \item{`read_macs2_broad`}
    \item{`read_seacr`}
    \item{`find_consensus`: find consensus peaks between two peak sets with adjust score, summit (MACS2), or maximun signals (SEACR)} 
    \item{`regionhit_per_sample_mat`}

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