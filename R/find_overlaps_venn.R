.get_hits_vectors <- function(x, y, hits) {
  # need a unit test for this
  hits_strings <- paste0(queryHits(hits), '-', subjectHits(hits))  
  
  left_hit <- paste0('left-', 1:length(x))
  right_hit <- paste0('right-', 1:length(y))

  left_hit[queryHits(hits)] <- hits_strings
  right_hit[subjectHits(hits)] <- hits_strings
  return(list(left_hit = left_hit, right_hit = right_hit))

}

#' @example
#' file_query <- system.file('extdata',
#'                           'chr2_Rep1_H1_H3K4me3.stringent.bed',
#'                           package = 'peaklerrr')
#' file_subject <- system.file('extdata',
#'                             'chr2_Rep2_H1_H3K4me3.stringent.bed',
#'                             package = 'peaklerrr')
#' file_subject <- system.file('extdata',
#'                             'chr2_Rep1_H1_CTCF.stringent.bed',
#'                             package = 'peaklerrr')
#'
#' query <- read_seacr(file_query)
#' subject <- read_seacr(file_subject)
#' @importFrom ggVennDiagram
#' @export 
find_overlaps_venn <- function(x, y, 
                               label_x=NULL,
                               label_y=NULL,
                               maxgap = -1L,
                               minoverlap = 1L, ...) {

  if (is.null(label_x)) label_x <- 'query'
  if (is.null(label_y)) label_y <- 'subject'
  
  hits <- findOverlaps(x, y, 
                       maxgap = maxgap, 
                       minoverlap = minoverlap, 
                       type = "any",
                       select = "all",
                       ignore.strand = TRUE)
  
  hits_vectors <- .get_hits_vectors(x, y, hits)
  names(hits_vectors) <- c(label_x, label_y)
  
  ggVennDiagram::ggVennDiagram(hits_vectors,
                               stroke_size = 0.5,
                               edge_lty = "blank",
                               edge_size = 0.5,
                               label = 'both', ...) +
    scale_fill_gradient(high = "#046C9A", low = "#ABDDDE" ) +
    theme(legend.position = 'none')
}