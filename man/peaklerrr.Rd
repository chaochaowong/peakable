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
    \item{`summit_from_seacr`}
    \item{`regionhit_per_sample_mat`}

  }
}

\section{Author}{
  \itemize{
    \item{Chao-Jen Wong (chaojen.wong@seattlechildrens.org)}
  }
}

\examples{
  \dontrun{
    data(expression_tpm)
    data(model)
    data(crispr_gene_dependency)
    data(crispr_gene_effect)
  }
}

\keyword{package}