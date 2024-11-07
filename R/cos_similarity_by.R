#' cos_similarity_by
#'
#' Given a sample sheet in data.frame and a list of peaks GRanges, find the
#' cos_similarity between peak sets, mostly replicates, grouped by
#' parameter sim_group_by'.
#'
#' @param sample_df A \code{data.frame} obtaining sample information
#' @param peaks_grl A list of peak ranges as \code{GRanges} instances
#' @param sim_group_by A vector of character to group the list of peak reanges
#'
#' @return A data.frame
#' @examples
#' add it later
#' @importFrom stringr str_replace_all
#' @importFrom purrr map_dfr
#' @export
cos_similarity_by <- function(sample_df,
                              peaks_grl,
                              sim_group_by,
                              return_data = TURE) {

  # validate sim_group_by
  is_valid <- all(sim_group_by %in% colnames(sample_df))
  invalid_cols <- sim_group_by[!sim_group_by %in% colnames(sample_df)]
  if (!is_valid)
    stop('Invalid columns in replicates_group_by:', invalid_cols)


  cos_sim <- .construct_cos_sim_grl(peaks_grl,
                                    sample_df,
                                    sim_group_by)

}

.paired_sim_df <- function(combn_pair, x, peaks_grl) {
  if (nrow(x) == 1) {
    sim_df <- data.frame(name = x$sample_id,
                         cos_sim =NA)
  }

  if (nrow(x) > 1 & all(x$number_of_peaks > 0)) {
    sim <-
      peakable::cos_similarity(peaks_grl[[x$sample_id[combn_pair[1]]]],
                               peaks_grl[[x$sample_id[combn_pair[2]]]])
    sim_df <- data.frame(name = paste(x$sample_id, collapse = '-vs-'),
                         cos_sim = sim)
  } else {
    sim_df <- data.frame(name = x$sample_id[1],
                         cos_sim =NA)
  }
  # need to carry the annotation
  return(sim_df)
}

.construct_cos_sim_grl <- function(peaks_grl,
                                   sample_df,
                                   sim_group_by) {
  sim_df <- sample_df %>%
    dplyr::group_split(across(all_of(sim_group_by))) %>%
    map(function(x) {
      m <- nrow(x)
      combn_pairs <- combn(c(1:m), m=m, simplify=FALSE)

      # get cos similarity between paired-samples in the group
      sim_df_by_group <- lapply(combn_pairs,
                                .paired_sim_df,
                                x, peaks_grl)
      # convert list of data.frame to a data.frame
      do.call(rbind, sim_df_by_group)
    })

  # set the name for sim_df, a list of data.frame
  keys <- sample_df %>%
    dplyr::group_by(across(all_of(sim_group_by))) %>%
    dplyr::group_keys()
  # set a unique name for each list element
  list_df <- keys %>%
    dplyr::mutate(group_id = apply(., 1, paste, collapse = "_")) %>%
    dplyr::mutate(group_id = stringr::str_replace_all(group_id,
  # set name for sim_df                                                     '_NA', ''))
  sim_df <- setNames(keys, list_df$group_id)
  # combine the list of data.frame with additional column indicating
  # how replicated are grouped together
  sim_df <- bind_rows(sim_df, .id='sim_group_by')
}
