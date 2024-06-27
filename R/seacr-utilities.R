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
#' @importFrom GenomeInfoDb dropSeqlevels keepStandardChromosomes seqlevelsStyle
#' @importFrom plyranges as_granges
#' @export
read_seacr <- function(file, drop_chrM = FALSE,
                       keep_standard_chrom = FALSE,
                       species = NULL) {
  # sanity check
  stopifnot(length(file) == 1)
  stopifnot(file.exists(file)) 
  
  # check if file is empty?
  is_empty <- file.info(file)$size == 0
  if (is_empty) return(NULL)
  

    # define column names
  col_names <- c("chr", "start", "end", "AUC", "max.signal", 
                 "max.signal.region", "num")

  # If file empty, return an empty GRange()
  is_empty <- file.info(file)$size == 0
  if (is_empty) return(GRanges())
  
  tb <- read.delim(file, header=FALSE)
  
  if (nrow(tb) > 0) {
    names(tb) <- col_names[1:ncol(tb)]
    tb$strand <- '*'
    gr <- plyranges::as_granges(tb, seqnames = chr)
    
    if (keep_standard_chrom & length(gr) > 0)
      gr <- GenomeInfoDb::keepStandardChromosomes(gr, 
                                                  pruning.mode='coarse',
                                                  species=species)
    
    if (drop_chrM & length(gr) > 0) {
      if (GenomeInfoDb::seqlevelsStyle(gr) == 'UCSC') value = 'chrM'
      if (GenomeInfoDb::seqlevelsStyle(gr) == 'NCBI') value = 'M'
      # if value (chrM or M) exists
      if (value %in% seqlevels(gr)) 
          gr <- GenomeInfoDb::dropSeqlevels(gr, value=value,
                                            pruning.mode='coarse')
    }
    
    return(gr)
  }
  
}

write_seacr <- function(x, file) {
  extra_cols <- c('AUC', 'max.signal',  'max.signal.region')
  valid_seacr <- all(extra_cols %in% names(mcols(x)))
  if (!valid_seacr) {
    stop(paste("For a valid SEACR peak there must be columns called:", 
               paste(extra_cols, collapse = ","), "in x."), call. = FALSE)
  }
  seacr_col_order <- c("seqnames", "start", "end", 
                    "AUC", "max.signal",  "max.signal.region")
  seacr_df <- as.data.frame(x)[, seacr_col_order]
  utils::write.table(seacr_df, file,
                     sep = "\t", row.names = FALSE,
                     col.names=FALSE,
                     na = ".", 
                     quote = FALSE)
}

#' Exact SEACR peak summit and convert to a GRanges object
#' @gr: A GRanges representing a peaks yield by SEACR, and it should include the 'max.signal.region' in the metadata columns
#' @summit_wid: NULL (default) or positive integer indicating the width of summit. Default to the width of 'max.signal.region'
#' 
#' @return a GRanges object with columns name and score (max.signal)
#' @rdname extract_summit
#' @importFrom plyranges mutate
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
  
  # santity check: mcols must max.signal.region
  flag <- all(c('AUC', 'max.signal', 'max.signal.region') %in% 
                names(mcols(gr)))
  if (flag %in% names(mcols(gr)))
    stop('The metadata columns must contain SEACR-specific columns AUC, max.signal, max.signal.region')
  
  # validate summit_wid
  if (!is.null(summit_wid) & !is.integer(summit_wid)) 
    summit_wid <- as.integer(summit_wid) 
    
  if (!is.null(summit_wid)) {
    if (summit_wid <= 0L) {
      message('Summit width must be greater than 0L. Default to NULL')
      summit_wid <- NULL
    }
  }
  
  summit <- GRanges(gr$max.signal.region) 
  summit <- plyranges::mutate(summit, 
                              name = paste0('peakname_', 1:length(gr)),
                              AUC = gr$AUC,
                              max.signal = gr$max.signal,
                              itemRgb ='#0000FF')
  
  # Ensure the width of summit conformed with the input setting
  if (!is.null(summit_wid) & is.integer(summit_wid)) {
    summit <- 
      plyranges::mutate(anchor_center(summit), width = summit_wid)
  }
  
  return(summit)
}