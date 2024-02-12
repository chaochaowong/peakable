#' utilities.R
#' TODO: all the utilities here need to be tested by unti test
#' 


#' @mat matrix (n-by-m)
#' @norm_factor a vector of length m
.normalize_columns <- function(mat, norm_factor) {
  #' multiple the m-th column of a N-by-M matrix by
  #' by the m-th element of norm_factor (1-by-M)
  stopifnot(ncol(mat) == length(norm_factor))
  t(t(mat) * norm_factor)
}

#' .getPCA()
#' @param mat, sample_info, n_pcs
#' @retrun data.frame with sample_id (if sample_info is provided),
#' PC1, PC2, ... PCn (depending on `n_pcs`), and columns
#' appended from sample_info
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join

.getPCA <- function(mat, sample_info=NULL, n_pcs=2) {
  pca <- prcomp(mat, scale. = TRUE)
  pcs <- as.data.frame(pca$rotation[, 1:n_pcs]) 
  
  if (!is.null(sample_info) & 'sample_id' %in% names(sample_info)) {
    pcs <- pcs %>%
      tibble::rownames_to_column(var='sample_id') %>%
      dplyr::left_join(sample_info, by='sample_id')
  }
  return(pcs)
}