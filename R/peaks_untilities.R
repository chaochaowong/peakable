consolidate_peaks <- function(grl) {
  # consolidate peaks from a list of GRanges
  # keep standard chromosomes
  # keep clean: no mcols and only standard chromosome
  gr <- unlist(as(grl, "GRangesList"))
  gr <- keepStandardChromosomes(gr, pruning.mode='coarse')
  mcols(gr) <- NULL
  GenomicRanges::reduce(gr)
}



