#' Import MACS2 narrowPeak and broadPeak BED file and format
#' ranges to a object
#' @file: a path to a file or a connection
#' @...: parameters pass to rtracklayer::import.bed()
#' 
#' @return a GRanges object
#' @rdname read_macs2
#' @importFrom rtracklayer import
#' @examples
#' narrow_file <- system.file('extdata', 
#'                            'chr2_Rep1_H1_CTCF_peaks.narrowPeak', 
#'                            package='peaklerrr')
#' gr <- read_macs2_narrow(narrow_file)
#' gr

# broadPeaks
#' broad_file <- system.file('extdata', 
#'                           'chr2_Rep1_H1_H3K27me3_peaks.broadPeak',
#'                            package='peaklerrr')
#' gr <- read_macs2_broad(broad_file)
#' gr
#' 
#' @export
read_macs2_narrow <- function(file, ...) {
  # sanity check
  stopifnot(length(file) == 1)
  stopifnot(file.exists(file)) 
  
  # define extra columns
  extraCols_narrowPeak <- c(signalValue = "numeric", 
                            pValue = "numeric",
                            qValue = "numeric", 
                            peak = "integer")
  
  rtracklayer::import(file, format = 'BED', 
                      extraCols = extraCols_narrowPeak, ...) 
}

#' @export
read_macs2_broad <- function(file, ...) {
  # sanity check
  stopifnot(length(file) == 1)
  stopifnot(file.exists(file)) 
  
  # define extra columns
  extraCols_broadPeak <- c(signalValue = "numeric", 
                           pValue = "numeric",
                           qValue = "numeric")
  
  rtracklayer::import(file, format = 'BED',
                      extraCols = extraCols_broadPeak, ...) 
}
