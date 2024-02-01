#' Import SEACR peak BED file and format
#' ranges to a object
#' @file: a path to a SEACR peak BED file or a connection
#' 
#' @return a GRanges object
#' @rdname read_seacr
#' @examples
#' seacr_file <- system.file('extdata',
#'                           'chr2_Rep1_H1_CTCF.stringent.bed',
#'                            package='peaklerrr')
#' gr <- read_seacr(seacr_file)
#' gr
#'
#' @export
read_seacr <- function(file) {
  # sanity check
  stopifnot(length(file) == 1)
  stopifnot(file.exists(file)) 
  
  # define column names
  col_names <- c("chr", "start", "end", "AUC", "max.signal", 
                 "max.signal.region", "num")

  tb <- read.delim(file)
  
  if (length(tb) > 0) {
    names(tb) <- col_names[1:ncol(tb)]
    tb$strand <- '*'
    plyranges::as_granges(tb, seqnames = chr)
  }
  
}

#' Exact SEACR peak summit and convert to a GRanges object
#' @gr: A GRanges representing a peaks yield by SEACR, and it should include the 'max.signal.region' in the metadata columns
#' @summit_wid: NULL (default) or positive integer indicating the width of summit. Default to the width of 'max.signal.region'
#' 
#' @return a GRanges object with columns name and score (max.signal)
#' @rdname extract_summit
#' @examples
#' seacr_file <- system.file('extdata',
#'                           'chr2_Rep1_H1_CTCF.stringent.bed',
#'                            package='peaklerrr')
#' gr <- read_seacr(seacr_file)
#' summit <- extract_summit_seacr(gr, summit_wid=1L)
#'
#' @export
extract_summit_seacr <- function(gr, summit_wid = NULL) {
  # convert to SEACR peak max.signal.region column to GRanges
  
  
  if (!is.null(summit_wid) & !is.integer(summit_wid)) 
    summit_wid <- as.integer(summit_wid) 
  if (!is.null(summit_wid) & summit_wid <= 0L) {
    message('Summit width must be greater than 0L. Default to NULL')
    summit_wid <- NULL
  }
  
  if ('max.signal.region' %in% names(mcols(gr)))
    stop('The metadata max.signal.region column does not exist.')
  
  summit <- GRanges(gr$max.signal.region) %>%
    plyranges::mutate(name = paste0('peakname_', 1:length(gr)),
                      score = gr$max.signal, 
                      itemRgb ='#0000FF')
  
  if (!is.null(summit_wid) & is.integer(summit_wid)) {
    summit <- 
      plyranges::mutate(anchor_center(summit), width = summit_wid)
  }
  
  return(summit)
}
