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

  if (keep_standard_chrom & length(gr) > 0)
    gr <- GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode='coarse',
                                                species = species)

  #if (drop_chrM & length(gr) > 0 ) {
  #  # this is only for human
  #  if (GenomeInfoDb::seqlevelsStyle(gr) == 'UCSC') value = 'chrM'
  #  if (GenomeInfoDb::seqlevelsStyle(gr) == 'NCBI') value = 'M'
  #  if (value %in% seqlevels(gr))
  #    gr <- GenomeInfoDb::dropSeqlevels(gr, value=value,
  #                                      pruning.mode='coarse')
  #}

  return(gr)
}

#' read_macs2_broad
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
  # check if the bed file has signalValue, pValue, qValue and not peak

  # what's the species? must assign species or genome
  gr <- rtracklayer::import(file, format = 'BED',
                      extraCols = extraCols_broadPeak, ...)

  # verify the sequence levels match the species

  if (keep_standard_chrom & length(gr) > 0)
    gr <- GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode='coarse',
                                                species = species)

  if (drop_chrM & length(gr) > 0) {

    if (GenomeInfoDb::seqlevelsStyle(gr) == 'UCSC') value = 'chrM'
    if (GenomeInfoDb::seqlevelsStyle(gr) == 'NCBI') value = 'M'
    if (value %in% seqlevels(gr))
      gr <- GenomeInfoDb::dropSeqlevels(gr, value=value,
                                        pruning.mode='coarse')
  }

  return(gr)
}

#' extract_summit_macs2
#'
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
    stop('gr must contain MACS2 narrowPeak specific column "peak".')

  if (is.null(summit_wid)) summit_wid <- 1L
  if (!is.integer(summit_wid)) summit_wid <- as.integer(summit_wid)
  if (summit_wid <= 0L)
    stop('Summit width must be greater than 0L.')

  # assuming the peak is not strand-sensitive
  summit <- gr %>%
    plyranges::mutate(start = start + peak, width=1L)
  summit <- plyranges::mutate(anchor_center(summit), width=summit_wid)
  mcols(summit) <- mcols(summit)[, names(mcols(summit)) != 'peak']

  return(summit)
}


#' write_broadpeaks
#'
write_broadpeaks <- function(x, file) {
  extra_cols <- c("signalValue", "pValue", "qValue")
  valid_broad <- all(extra_cols %in% names(mcols(x)))

  if (!valid_broad) {
    stop(paste("For a valid MACS2 broad peak there must be columns called:",
               paste(extra_cols, collapse = ","), "in x."), call. = FALSE)
  }

  broad_col_order <- c("seqnames", "start", "end", "name", "score",
                       "strand", extra_cols)
  broad_df <- as.data.frame(x)[, broad_col_order]

  utils::write.table(broad_df, file,
                     sep = "\t", row.names = FALSE,
                     col.names=FALSE,
                     na = ".",
                     quote = FALSE)
}
