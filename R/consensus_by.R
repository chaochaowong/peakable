#' @importFrom stringr str_detect
#'
consensus_by <- function(sample_df, peaks_grl,
                         consensus_group_by,
                         peak_caller = c('macs2', 'searc')
                        ) {
  # wrapper function of find_consensus_macs2 and find_consensus_searc
  # consensus_by = c('cell_line', 'antibody')
  consensus <- vector('list', 2)

  # validate consensus_group_by
  is_valid <- all(consensus_group_by %in% colnames(sample_df))
  invalid_cols <- consensus_group_by[!consensus_group_by %in% colnames(sample_df)]
  if (!is_valid)
        stop('Invalid columns in consensus_group_by:', invalid_cols)

  # validate peak_caller
  if (!is.null(peak_caller)) {
    peak_caller <- tolower(peak_caller[1])

    if (!str_detect(peak_caller, 'seacr|macs2'))
      peak_call <- NULL
  }

  # define overlap functions
  if (str_detect(peak_caller, 'seacr'))
    overlap_call_func <- peakable::find_consensus_seacr
  if (str_detect(peak_caller, 'macs2'))
    overlap_call_func <- peakable::find_consensus_macs2
  if (is.null(peak_caller))
    overlap_call_func <- plyranges::find_overlaps

  # construct consensus$grl using consensus_group_by
  consensus <- .construct_consensus_grl(peaks_grl,
                                        sample_df,
                                        consensus_group_by,
                                        overlap_call_func)
  # consensus$df includes columns:
  #  - sample_id
  #  - number_of_peaks
  #  - append the consensus_group_b columns
  return(consensus)
}

.construct_consensus_grl <- function(peaks_grl, sample_df,
                                     consensus_group_by,
                                     overlap_call_func) {
  # can only handle two replicates
  grl <- sample_df %>%
    dplyr::group_split(across(all_of(consensus_group_by))) %>%
    map(function(df) {
      if (nrow(df) == 1) { peaks_grl[[df$sample_id[1]]] }
      if (nrow(df) > 1) {
        x <- peaks_grl[[df$sample_id[1]]]
        y <- peaks_grl[[df$sample_id[2]]]
        overlap_call_func(x, y)
      }
   })

  keys <- sample_df %>%
    dplyr::group_by(across(all_of(consensus_group_by))) %>%
    dplyr::group_keys()

  # Create a unique name for each list element
  list_df <- keys %>%
    dplyr::mutate(sample_id = apply(., 1, paste, collapse = "_"))


  # granges list
  grl <- setNames(grl, list_df$sample_id)

  grl_df <- data.frame(sample_id = names(grl),
                   number_of_peaks = sapply(grl, length)) %>%
    dplyr::left_join(list_df, by='sample_id')

  return(list(grl=grl, df=grl_df))
}

