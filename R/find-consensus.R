#' Find consensus of two SEACR peaks: a wrapping function of
#' plyranges::find_overlaps() and modify the metadata columns such as
#' AUC, max.signal and max.signal.regions
#' @x right object representing SEACR peaks
#' @y left object representing SEACR peaks
#' @minoverlap NULL
#'
#' @return a GRanges object
#' @rdname read_seacr
#' @examples
#' seacr_file_x <-
#'   system.file('extdata',
#'              'chr2_Rep1_H1_CTCF.stringent.bed',
#'               package='peakable')
#' seacr_file_y <-
#'   system.file('extdata',
#'              'chr2_Rep2_H1_CTCF.stringent.bed',
#'                package='peakable')
#'
#' x <- read_seacr(seacr_file_x)
#' y <- read_seacr(seacr_file_y)
#' consensus <- find_consensus_seacr(x, y, minoverlap = 40L)
#' consensus
#'
#' consensus <- find_consensus_seacr(x, y)
#' consensus
#' @importFrom plyranges find_overlaps
#' @export
find_consensus_seacr <- function(x, y, minoverlap = NULL) {

  # sanity check; X and Y must be GRanges
  stopifnot(length(x) >= 1)
  stopifnot(length(y) >= 1)

  if (is.null(minoverlap)) {
    min_wid <- min(min(width(x)), min(width(y)))
    minoverlap <- min_wid / 2
  }

  # use plyranges::find_overlap
  consensus <-
    plyranges::find_overlaps(x, y, minoverlap = minoverlap)

  # update mcols: max.signal.region, max.signal, AUC
  consensus_mcols <- mcols(consensus)
  keep_x <- consensus_mcols$max.signal.x > consensus_mcols$max.signal.y

  consensus$AUC <- ifelse(keep_x, consensus_mcols$AUC.x,
                          consensus_mcols$AUC.y)
  consensus$max.signal <- ifelse(keep_x, consensus_mcols$max.signal.x,
                                 consensus_mcols$max.signal.y)
  consensus$max.signal.region <- ifelse(keep_x,
                                        consensus_mcols$max.signal.region.x,
                                        consensus_mcols$max.signal.region.y)

  # drop off the mcols of x and y and keep unique
  mcols(consensus) <- mcols(consensus)[, c('AUC', 'max.signal',
                                           'max.signal.region'),
                                       drop = FALSE]
  BiocGenerics::unique(consensus)
}

#' find_consensus_macs2
#'
#' Find consensus of two MACS2 peaks: a wrapping function of
#' plyranges::find_overlaps() and modify the metadata columns such as
#'
#' @param x right object representing MACS2 narrow peaks
#' @param y left object representing MACS2 narrow peaks
#' @param minoverlap NULL
#'
#' @return a GRanges object
#' @rdname find_consensus_macs2
#' @examples
#' macs2_file_x <-
#'   system.file('extdata',
#'              'chr2_Rep1_H1_CTCF_peaks.narrowPeak',
#'               package='peakable')
#' macs2_file_y <-
#'   system.file('extdata',
#'              'chr2_Rep2_H1_CTCF_peaks.narrowPeak',
#'                package='peakable')
#'
#' x <- read_macs2_narrow(macs2_file_x)
#' y <- read_macs2_narrow(macs2_file_y)
#' consensus <- find_consensus_macs2(x, y, minoverlap = 40L)
#' consensus
#'
#' consensus <- find_consensus_macs2(x, y)
#' consensus
#' @importFrom plyranges find_overlaps
#' @importFrom stringr str_detect
#' @export
find_consensus_macs2 <- function(x, y, minoverlap = NULL) {
  # sanity check; X and Y must be GRanges
  stopifnot(length(x) >= 1)
  stopifnot(length(y) >= 1)

  if (is.null(minoverlap)) {
    min_wid <- min(min(width(x)), min(width(y)))
    minoverlap <- min_wid / 2
  }

  # use plyranges::find_overlap
  consensus <-
    plyranges::find_overlaps(x, y, minoverlap = minoverlap)

  # only keep mcols of x and unique only; letter will change to whichever have
  # better qValue
  consensus_mcols <- mcols(consensus)
  consensus$name <- consensus_mcols$name.x
  consensus$score <- consensus_mcols$score.x
  consensus$signalValue <- consensus_mcols$signalValue.x
  consensus$pValue <- consensus_mcols$pValue.x
  consensus$qValue <- consensus_mcols$qValue.x

  keep_cols <- c('name', 'score', 'signalValue', 'pValue', 'qValue')
  if (any(stringr::str_detect(names(mcols(consensus)), 'peak'))) {
    consensus$peak <- consensus_mcols$peak.x
    keep_cols <- c(keep_cols, 'peak')
  }

  mcols(consensus) <- mcols(consensus)[, keep_cols, drop = FALSE]
  return(BiocGenerics::unique(consensus))
}
