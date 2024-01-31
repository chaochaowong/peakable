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



summit_from_seacr <- function(seacr_gr, summit_wid = NULL) {
  # convert to seacr max.signal.region column to GRanges
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
