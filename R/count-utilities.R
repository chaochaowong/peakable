#' peak_read_count() count reads hitting the features (GRanges). If
#' spike_in_norm is true, normalize the count matrix by spike_in_factor.
#' @usage peak_read_count(features, sample_df, spike_in_norm = FALSE)
#' @features GRanges object of features
#' @sample_df data.frame object of sample information
#' @spike_in_norm logical indicating whether to normalized by spike_in_factor
#' @retrun RangeSummarizedExperiment
#' @return a RangesSummarizedExperiment object
#' @importFrom Rsamtools BamFileList
#' @importFrom GenomicAlignments summarizeOverlaps
#' @example 
#' 
peak_read_count <- function(features, sample_df, 
                            spike_in_norm = FALSE) {
  # - sample_df is a data.frame must have sample_id, bam_file, aligned_paired, 
  #   and spike_in_factor (spike_nrom=TRUE) columns
  # - count reads hitting the peakset then normlized by spike_in_norm (if spike_norm=TRUE)
  stopifnot('sample_id' %in% names(sample_df))
  stopifnot('bam_file' %in% names(sample_df))
  stopifnot(is(features, 'GRanges'))
  
  # sanity check: file existence; check sample_df columns
  bam_files <- Rsamtools::BamFileList(sample_df$bam_file)
  
  se <- 
    GenomicAlignments::summarizeOverlaps(features = features,
                                         reads = bam_files,
                                         ignore.strand = TRUE,
                                         singleEnd = FALSE,
                                         mode = "Union",
                                         fragments=FALSE,
                                         inter.feature=FALSE)
  # tidy rownames and column data of se
  rownames(se) <- paste0('peakname_', 1:length(se))
  colData(se) <- DataFrame(sample_df)
  colnames(se) <- sample_df$sample_id
  
  se$reads <- colSums(assays(se)[['counts']])
  
  #' include spike_in_norm assay
  if (spike_in_norm) {
    if ('spike_in_factor' %in% names(colData(se))) {
      assays(se)[['spike_in_norm']] <- 
        .normalize_columns(mat = assays(se[['count']]),
                           norm_factor = se$spike_in_factor)
      # normalized "reads" by spike-in factor
      se$reads_spikein_norm <- se$reads * se$spike_in_factor
    }
  }
  
  #' if read_paired (lib size) is available,
  #' then give FRiP score and lib size factor for
  #' linear normalization based on library size (read_paired)
  if ('read_paired' %in% names(colData(se))) {
    se <- .estimate_FRiP_score(se) # Fraction of reads overlap with features
    se <- .estimate_lib_size_factor(se) # can be used for normalization
  }
  
  return(se)
  
}

