#' peakle_flow.R
#' @example 
#' dir_header <- '///Volumes' # via MacBook
#' result_dir <- file.path(dir_header, 'Active',
#'                         'lawlor_e/Shireen/CnT_Results',
#'                         'A673_anthracyclines_240130')
#' sample_df <- read_csv(file.path(result_dir, 'data', 'nf-sample-sheet.csv'))                         
#' 
#' peak_bed_dir <- file.path(result_dir, 'peaks_calls', 'seacr_callpeak')
#' peak_bed_pattern <- '\\_threshold0.01_non.stringent.bed$'
#' seacr <- 
#'   peaklerrr:::peakle_flow(sample_df, 
#'                           result_dir, result_dir,
#'                           peak_caller = 'SEACR-thres1p',
#'                           peak_bed_dir = peak_bed_dir,
#'                           peak_bed_pattern = peak_bed_pattern)
#'
#'  seacr_hit_mat <- peaklerrr::consolidated_peak_hits(seacr$grl)
#' seacr_hit_pca <- 
#'   peaklerrr:::.getPCA(seacr_hit_mat, 
#'                       sample_info=seacr$df, n_pcs=2)
#' ggplot(seacr_hit_pca, aes(x=PC1, y=PC2, 
#'                          color=antibody, shape=treatment)) +
#'   geom_point() + theme_minimal() +
#'   labs(title='PCA: SEACR peak-hits matrix')
#'   
.set_bam_params <- function(result_dir, bam_pattern) {
  #' ignore bam_dir and bam_pattern
  bam_dir <- file.path(result_dir, 'samtools_sort')

  if (is.null(bam_pattern)) {
    bam_pattern <- '\\.markedDup.filter.sort.bam$'
  } 
  
  #' sanity checK does the bam file exists?
  return(list(bam_dir = bam_dir, bam_pattern = bam_pattern))
}

.set_samtools_params <- function(result_dir) {
  #' ignore bam_dir and bam_pattern
  stats_dir <- file.path(result_dir, 'samtools_stats')
  stats_pattern <- '\\.markedDup.stats$'
  #' sanity checK does stats_dir exists?
  return(list(stats_dir = stats_dir, stats_pattern = stats_pattern))
}

.get_list_files_and_sample_id <- function(path, pattern) {
  data.frame(
    file = list.files(path, pattern = pattern,
                      full.names=TRUE)) %>%
    dplyr::mutate(sample_id = str_replace(basename(file),
                                          pattern, ''))
}

.samtools_stats <- function(path, pattern) {
  if (!file.exists(path)) stop('samtools_stats path does not exist.')
  
  # get duplication rate
  
  stats_df <- .get_list_files_and_sample_id(path, pattern)
  
  # parse the log
  stats_log <- purrr::map_dfr(stats_df$file, function(x) {
    ln <- readr::read_delim(x, delim="\t", skip=16, 
                            n_max=3, col_names=FALSE)[, -1]
    data.frame(sample_id = str_replace(basename(x), '.markedDup.stats', ''),
               read_paired = ln[2, 'X3', drop=TRUE],
               read_duplicated = ln[3, 'X3', drop=TRUE]) }) %>%
    dplyr::mutate(dup_rate = read_duplicated / read_paired)
}

peakle_flow <- function(sample_df, # must be from nf_sample_sheet
                        result_dir,
                        peak_caller, 
                        peak_bed_dir, 
                        peak_bed_pattern,
                        bam_pattern = NULL) {
  library(tidyverse)
  require(purrr)
  require(BiocParallel)
  
  #' setting parameters:
  #' if result_dir is given, then use default bam directory
  if (!file.exists(result_dir))
    stop(result_dir, ' does not exist.')
  
  if (file.exists(result_dir)) {
      
      bam_params <- .set_bam_params(result_dir, bam_pattern)
      bam_dir <- bam_params$bam_dir
      bam_pattern <- bam_params$bam_pattern
      samtools_params <- .set_samtools_params(result_dir)
      stats_dir <- samtools_params$stats_dir
      stats_pattern <- samtools_params$stats_pattern
  }

  # check the input sample_df is valid
  # check if stats_dir exists, if not, skip samtools stats
  
  # 1.a) bam files; sample_id is the basename with bam_patter replaced
  bam_df <- .get_list_files_and_sample_id(bam_dir, bam_pattern) %>%
    dplyr::rename(bam_file = 'file')
  
  # check if stats_dir exists
  if (file.exists(stats_dir) ) {
     stats_df <- .samtools_stats(stats_dir, stats_pattern)
  } else {
    stats_df <- NULL
  }
     
  
  # 2) append to sample_df; remove irrelavent columns
  sample_df <- sample_df %>%
    dplyr::left_join(bam_df, by='sample_id') 
  
  if (!is.null(stats_df)) {
    sample_df <- sample_df %>%
      dplyr::left_join(stats_df, by='sample_id')
  }
  
  # 3) get peak bed files
  peak_df <- data.frame(
    bed_file = list.files(peak_bed_dir, pattern=peak_bed_pattern,
                          full.name=TRUE)) %>%
    dplyr::mutate(peakcall_id = str_replace(basename(bed_file),
                                          peak_bed_pattern, ''),
                  peak_caller = peak_caller) %>%
    # if relaxed or stringent compared with IgG, then tidy sample_id
    dplyr::mutate(sample_id =  str_split(peakcall_id, '_vs_', 
                                         simplify=TRUE)[, 1]) %>%
    dplyr::right_join(sample_df, by='sample_id') %>%
    dplyr::filter(!is.na(bed_file)) # some IgG files do not have bed files
  
  # additional column for the IgG file
  if (str_detect(peak_caller, 'relaxed|stringent')) {
    peak_df <- peak_df %>%
      dplyr::mutate(IgG = str_split(peakcall_id, '_vs_', 
                                    simplify=TRUE)[, 2])
  }
  
  # tidyup and remove peakcall_id
  peak_df <- peak_df %>% dplyr::select(-peakcall_id)
  
  # 4) sanity check
  if (nrow(peak_df) < 1)
    stop('Cannot not fine the peak bed files for these samples.')
  # 4) define the peak ranges

  if (str_detect(peak_caller, 'SEACR')) 
    peak_call_func <- peaklerrr::read_seacr
  if (str_detect(peak_caller, 'narrow'))
    peak_call_func <- peaklerrr::read_macs2_narrow
  if (str_detect(peak_caller, 'broad'))
    peak_call_func <- peaklerrr::read_macs2_broad
  
  message('Peak callers: ', peak_caller)
  message('Get peak ranges ...')
  
  peak_grl <- bplapply(peak_df$bed_file, 
                       peak_call_func,
                       drop_chrM = TRUE,
                       keep_standard_chrom = TRUE,
                       species = 'Homo_sapiens')
  names(peak_grl) <- peak_df$sample_id
  
  peak_df$number_of_peaks <- sapply(peak_grl, length)
  
  return(list(df = peak_df, grl = peak_grl))
}