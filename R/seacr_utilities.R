summit_from_seacr <- function(seacr_gr, summit_wid = NULL) {
  # convert to seacr max.signal.region column to GRanges
  summit <- GRanges(gr$max.signal.region) %>%
    plyranges::mutate(name = paste0('peakname_', 1:length(sub)),
                      score = gr$max.signal, 
                      itemRgb ='#0000FF')
  if (!is.null(summit_wid) & is.integer(summit_wid)) {
    summit <- 
      plyranges::mutate(anchor_center(summit), width = summit_wid)
  }
  return(summit)
}

read_seacr <- function(file) {
  # sanity check
  stopifnot(length(file) == 1)
  stopifnot(file.exists(file)) 
  
  # define column names
  col_names <- c("chr", "start", "end", "AUC", "max.signal", 
                 "max.signal.region", "num")
  
  tb <- read.delim(file, col.names = FALSE)
  
  if (length(tb) > 0) {
    names(tb) <- col_names[1:ncol(tb)]
    tb$strand <- '*'
    gr <- plyranges::as_granges(tb, seqnames = chr)
  }
  
}