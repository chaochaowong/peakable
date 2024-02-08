#' Import MACS2 narrowPeak and broadPeak BED file and format
#' ranges to a object
#' @file: a path to a file or a connection
#' @...: parameters pass to rtracklayer::import.bed()
#' 
#' @return a GRanges object
#' @rdname read_macs2
#' @importFrom rtracklayer import
#' @importFrom GenomeInfoDb dropSeqlevels keepStandardChromosomes seqlevelsStyle
#' @examples
#' narrow_file <- system.file('extdata', 
#'                            'chr2_Rep1_H1_CTCF_peaks.narrowPeak', 
#'                            package='peaklerrr')
#' gr <- read_macs2_narrow(narrow_file)
#' gr
#' @export
read_macs2_narrow <- function(file, drop_chrM = FALSE, 
                              keep_standard_chrom = FALSE,
                              species = NULL, ...) {
  # sanity check
  stopifnot(length(file) == 1)
  stopifnot(file.exists(file)) 
  
  # define extra columns
  extraCols_narrowPeak <- c(signalValue = "numeric", 
                            pValue = "numeric",
                            qValue = "numeric", 
                            peak = "integer")
  
  gr <- rtracklayer::import(file, format = 'BED', 
                      extraCols = extraCols_narrowPeak, ...) 
  
  if (keep_standard_chrom)
    gr <- GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode='coarse',
                                                species = species)
  
  if (drop_chrM) {
    if (GenomeInfoDb::seqlevelsStyle(gr) == 'UCSC') value = 'chrM'
    if (GenomeInfoDb::seqlevelsStyle(gr) == 'NCBI') value = 'M'
    if (value %in% seqlevels(gr)) 
      gr <- GenomeInfoDb::dropSeqlevels(gr, value=value,
                                        pruning.mode='coarse')
  }
  return(gr)
}

#' Import MACS2 broadPeak and broadPeak BED file and format
#' ranges to a object
#' @file: a path to a file or a connection
#' @...: parameters pass to rtracklayer::import.bed()
#' 
#' @return a GRanges object
#' @rdname read_macs2
#' @importFrom rtracklayer import
#' @importFrom GenomeInfoDb dropSeqlevels keepStandardChromosomes seqlevelsStyle
#' @example
#' broadPeaks
#' broad_file <- system.file('extdata', 
#'                           'chr2_Rep1_H1_H3K27me3_peaks.broadPeak',
#'                            package='peaklerrr')
#' gr <- read_macs2_broad(broad_file)
#' gr
#' 
#' @export
read_macs2_broad <- function(file, drop_chrM = FALSE, 
                             keep_standard_chrom = FALSE,
                             species = NULL, ...) {
  # sanity check
  stopifnot(length(file) == 1)
  stopifnot(file.exists(file)) 
  
  # define extra columns
  extraCols_broadPeak <- c(signalValue = "numeric", 
                           pValue = "numeric",
                           qValue = "numeric")
  
  gr <- rtracklayer::import(file, format = 'BED',
                      extraCols = extraCols_broadPeak, ...) 
  
  if (keep_standard_chrom)
    gr <- GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode='coarse',
                                                species = species)
  
  if (drop_chrM) {
    if (GenomeInfoDb::seqlevelsStyle(gr) == 'UCSC') value = 'chrM'
    if (GenomeInfoDb::seqlevelsStyle(gr) == 'NCBI') value = 'M'
    if (value %in% seqlevels(gr)) 
      gr <- GenomeInfoDb::dropSeqlevels(gr, value=value,
                                        pruning.mode='coarse')
  }
  return(gr)
}

#' Exact MACS2 narrow peak summit and convert to a GRanges object
#' @gr: A GRanges of MACS2 narrow peak, and it should include the 'max.signal.region' in the metadata columns
#' @summit_wid: NULL (default) or positive integer indicating the width of summit. Default to the width of 'max.signal.region'
#' 
#' @return a GRanges object
#' @rdname extract_summit
#' @importFrom plyranges as_granges
#' @examples
#' narrow_file <- system.file('extdata', 
#'                            'chr2_Rep1_H1_CTCF_peaks.narrowPeak', 
#'                            package='peaklerrr')
#' gr <- read_macs2_narrow(narrow_file)
#' summit <- extrac_summit_macs2(gr, summit_wid=1L)
#' summit
#' @export
extract_summit_macs2 <- function(gr, summit_wid = NULL) {
  
  # check mcols are MACS2-specific
  if (!'peak' %in% names(mcols(gr)))
    stop('gr must contain MACS2-specific column "peak".')
  
  if (is.null(summit_wid)) summit_wid <- 1L
  if (!is.integer(summit_wid)) summit_wid <- as.integer(summit_wid) 
  if (summit_wid <= 0L)
    stop('Summit width must be greater than 0L.')
  
  
  summit <- data.frame(start=start(gr) + gr$peak,
                       width = summit_wid,
                       seqnames = seqnames(gr),
                       strand = '*')
  summit <- plyranges::as_granges(summit)
  mcols(summit) <- mcols(gr)[, names(mcols(gr)) != 'peak']
  summit
}
