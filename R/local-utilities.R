#' utilities.R
#' TODO: all the utilities here need to be tested by unti test
#' 

.normalize_columns <- function(mat, norm_factor) {
  #' multiple the m-th column of a N-by-M matrix by
  #' by the m-th element of norm_factor (1-by-M)
  stopifnot(ncol(mat) == length(norm_factor))
  t(t(mat) * norm_factor)
}

#' ====
#' .getPCA()
#' @param mat, sample_info, n_pcs
#' @retrun data.frame with sample_id (if sample_info is provided),
#' PC1, PC2, ... PCn (depending on `n_pcs`), and columns
#' appended from sample_info
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join

.getPCA <- function(mat, sample_info=NULL, n_pcs=2) {
  pca <- prcomp(t(mat))
  pcs <- as.data.frame(pca$x[, 1:n_pcs]) 
  
  if (!is.null(sample_info) & 'sample_id' %in% names(sample_info)) {
    pcs <- pcs %>%
      tibble::rownames_to_column(var='sample_id') %>%
      dplyr::left_join(sample_info, by='sample_id')
  }
  
  var_pcs <- pca$sdev^2
  # proportion of total variance explained
  prop_var <- var_pcs / sum(var_pcs)
  return(list(pcs = pcs, prop_var = prop_var[1:n_pcs] )
}

#' ====
#' @param se A RangedSummarizedExperiment or DESeqDataSet object
.estimate_FRiP_score <- function(se) {
  stopifnot(is(se, "RangedSummarizedExperiment"))
  se$FRiP <- se$reads / se$read_paired
  
  return(se)
}

#' ====
#' @param se A RangedSummarizedExperiment or DESeqDataSet object
.estimate_lib_size_factor <- function(se) {
  stopifnot(is(se, "RangedSummarizedExperiment"))
  
  if ('read_paired' %in% names(colData(se))) {
    lib_mean <- mean(se$read_paired)
    se$lib_size_factor <- se$read_paired / lib_mean
  }
  
  return(se)
}

#' ====
#' .peak_coverage_pearson() generates Pearson correlation heatmap for 
#' peak read coverage. The read coverage object should be a DESeqTransform
#' object. If it is a DESeqDataSet object, it will be transform into
#' a DESeqTransfrom object by DESeq2::rlog(). The .peak_coverage_pearson()
#' employs corrr::correlate() to get the Pearson coefficient and uses 
#' pheatmap::pheatmap() to generate the heatmap. The \code{...} argument
#' are the parameters to pass to pheatmap::pheatmap().
#' 
#' @param se either a DESeqTransform (preferred) or DESeqDataSet object
#' @param file_name file path where to save the figure
#' @param ... arguments to pass to pheatmap::pheatmap()
#' @importFrom corrr correlate
#' @importFrom dplyr mutate
#' @importFrom DESeq2 rlog
#' @importFrom pheatmap pheatmap
.peak_coverage_pearson <- function(se, file_name = NULL, ...) {
  require(corrr)
  require(pheatmap)
  require(dplyr)
  
  #' rlg must be DESeqTransform object
  stopifnot(is(se, 'DESeqTransform') | is(se, 'DESeqDataSet'))
  
  #' if not DESeqTransfrom, transform by "rlog"
  if (!is(se, 'DESeqTransform'))
    se <- DESeq2::rlog(se, blind = TRUE)
  
  se_cor <- assay(se) %>%
    corrr::correlate() %>%
    dplyr::mutate(across(everything(), ~replace_na(.x, 1)))
  
  if (!is.null(file_name))
    pheatmap(as.matrix(se_cor[, -1]), show_colnames=TRUE,
             angle_col=90, fontsize_col=8, width=8,
             scale = "none",
             display_numbers = TRUE,
             number_format = '%.2f',
             fontsize_number = 5,
             slient = TRUE,
             filename = file_name, ...)
  
  if (is.null(file_name))
    pheatmap(as.matrix(se_cor[, -1]), show_colnames=TRUE,
             angle_col=90, fontsize_col=8, width=8,
             scale = "none",
             display_numbers = TRUE,
             number_format = '%.2f',
             fontsize_number = 5)
  
}