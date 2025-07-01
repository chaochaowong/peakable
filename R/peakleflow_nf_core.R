#' peakle_flow.R
#' @example
#' /donotrun{
#' dir_header <- '///Volumes' # via MacBook
#' result_dir <- file.path(dir_header, 'Active',
#'                         'lawlor_e/Shireen/CnT_Results',
#'                         'A673_anthracyclines_240130')
#' sample_df <- read_csv(file.path(result_dir, 'data', 'nf-sample-sheet.csv'))
#'
#' peak_bed_dir <- file.path(result_dir, 'peaks_calls', 'seacr_callpeak')
#' peak_bed_pattern <- '\\_threshold0.01_non.stringent.bed$'
#' seacr <-
#'   peakable:::peakle_flow(sample_df,
#'                           result_dir, result_dir,
#'                           peak_caller = 'SEACR-thres1p',
#'                           peak_bed_dir = peak_bed_dir,
#'                           peak_bed_pattern = peak_bed_pattern)
#'
#'  seacr_hit_mat <- peakable::consolidated_peak_hits(seacr$grl)
#' seacr_hit_pca <-
#'   peaklable:::.getPCA(seacr_hit_mat,
#'                       sample_info=seacr$df, n_pcs=2)
#' ggplot(seacr_hit_pca, aes(x=PC1, y=PC2,
#'                          color=antibody, shape=treatment)) +
#'   geom_point() + theme_minimal() +
#'   labs(title='PCA: SEACR peak-hits matrix')
#'}
.set_bam_params_nf_core <- function(result_dir, bam_dir, bam_pattern) {
  #' ignore bam_dir and bam_pattern
  if (is.null(bam_dir))
    bam_dir <- file.path(result_dir, '02_alignment/bowtie2/target/markdup')

  if (is.null(bam_pattern)) {
    bam_pattern <- '\\.target.markdup.sorted.bam$'
  }

  #' sanity checK does the bam file exists?
  return(list(bam_dir = bam_dir, bam_pattern = bam_pattern))
}

.set_samtools_params_nf_core <- function(result_dir) {
  #' ignore bam_dir and bam_pattern
  stats_dir <- file.path(result_dir, '02_alignment/bowtie2/target/markdup')
  stats_pattern <- '\\.stats$'
  #' sanity checK does stats_dir exists?
  return(list(stats_dir = stats_dir, stats_pattern = stats_pattern))
}

.get_list_files_and_sample_id_nf_core <- function(path, pattern) {
  data.frame(
    file = list.files(path, pattern = pattern,
                      full.names=TRUE)) %>%
    dplyr::mutate(sample_id = str_replace(basename(file),
                                          pattern, ''))
}

.samtools_stats_nf_core <- function(path, pattern) {
  if (!file.exists(path)) stop('samtools_stats path does not exist.')

  # get duplication rate

  stats_df <- .get_list_files_and_sample_id_nf_core(path, pattern)

  # parse the log
  stats_log <- purrr::map_dfr(stats_df$file, function(x) {
    ln <- readr::read_delim(x, delim="\t", skip=16,
                            n_max=3, col_names=FALSE)[, -1]
    data.frame(sample_id = str_replace(basename(x), '.stats', ''),
               read_paired = ln[2, 'X3', drop=TRUE],
               read_duplicated = ln[3, 'X3', drop=TRUE]) }) %>%
    dplyr::mutate(dup_rate = read_duplicated / read_paired)
}

peakleflow_nf_core <- function(sample_df, # must be from nf_sample_sheet
                               result_dir,
                               peak_caller,
                               peak_bed_dir,
                               peak_bed_pattern,
                               bam_dir = NULL, 
                               bam_pattern = NULL, 
                               species = "Homo_sapiens",
                               remove_blacklist = FALSE,
                               blacklist_file = NULL,
                               drop_chrM = FALSE
                               mito = NULL) {
  require(dplyr)
  require(stringr)
  require(purrr)
  require(BiocParallel)

  #' setting parameters:
  #' if result_dir is given, then use default bam directory
  if (!file.exists(result_dir))
    stop(result_dir, ' does not exist.')

  if (!file.exists(peak_bed_dir))
    stop(peak_bed_dir, ' does not exist.')

  if (!is.null(bam_dir)) {
     stopifnot(file.exists(bam_dir))
  }

  # define bam_dir and patterns
  bam_params <- .set_bam_params_nf_core(result_dir, bam_dir, bam_pattern)

  # define stats_dir
  samtools_params <- .set_samtools_params_nf_core(result_dir)

  # assign parameters
  bam_dir <- bam_params$bam_dir
  bam_pattern <- bam_params$bam_pattern
  stats_dir <- samtools_params$stats_dir
  stats_pattern <- samtools_params$stats_pattern
  peak_caller <- tolower(peak_caller[1])

  # check the input sample_df is valid
  # check if stats_dir exists, if not, skip samtools stats

  # 1.a) bam files; sample_id is the basename with bam_patter replaced
  bam_df <- .get_list_files_and_sample_id_nf_core(bam_dir, bam_pattern) %>%
    dplyr::rename(bam_file = 'file')

  if (nrow(bam_df) == 0)
    stop('BAM file pattern does not exist: ', bam_pattern)
  # append bam_df to sample_df by 'sample_id'

  sample_df <- sample_df %>%
    dplyr::inner_join(bam_df, by='sample_id')

  if (nrow(sample_df) == 0)
    stop('The sample_id extracted from BAM files does not match sample_df$sample_id.')

  # check if stats_dir exists
  if (file.exists(stats_dir) ) {
    message('Sorting samtools stats ...')
    stats_df <- .samtools_stats_nf_core(stats_dir, stats_pattern)
  } else {
    message(stats_dir, ' does not exist. Set samtools stats to null')
    stats_df <- NULL
  }


  # 2) append to stats_df to sample_df
    if (!is.null(stats_df)) {
    sample_df <- sample_df %>%
      dplyr::left_join(stats_df, by='sample_id')
  }


  # 3) get peak bed files
  bed_files <- list.files(peak_bed_dir, pattern=peak_bed_pattern,
                          full.name=TRUE)
  if (identical(bed_files, character(0))) {
    msg <- sprintf('Cannot find matting bed pattern %s in %s',
                   peak_bed_pattern, peak_bed_dir)
    stop(msg)
  }

  # append all bed_files to sample_df and other information
  peak_df <- data.frame(bed_file = bed_files) %>%
    dplyr::mutate(peakcall_id = str_replace(basename(bed_file),
                                            peak_bed_pattern, ''),
                  peak_caller = peak_caller) %>%
    # if relaxed or stringent compared with IgG, then tidy sample_id
    dplyr::mutate(sample_id =  str_split(peakcall_id, '_vs_',
                                         simplify=TRUE)[, 1]) %>%
    dplyr::inner_join(sample_df, by='sample_id') %>%
    dplyr::filter(!is.na(bed_file)) # some IgG files do not have bed files

  if (nrow(peak_df) < 1)
    stop('Cannot not match the peak bed files and bam files for these samples.')

  # additional column for the IgG file: if caller is SEACR
  if (str_detect(peak_caller, 'seacr')) { # not detect threshold
    if (str_detect(peak_caller, 'relaxed|stringent')) {
       peak_df <- peak_df %>%
        dplyr::mutate(IgG = str_split(peakcall_id, '_vs_',
                                      simplify=TRUE)[, 2])
    }
  }

  #
  # check point
  #
  if (!identical(nrow(peak_df), nrow(sample_df))) {
    msg <- sprintf('The number of peak files does not match with that of the non-IgG samples: %s',
                   setdiff(sample_df$sample_id, peak_df$sample_id))
    message(msg)
  }

  # tidyup and remove peakcall_id
  peak_df <- peak_df %>% dplyr::select(-peakcall_id)

  # 4) define the peak ranges reader

  if (str_detect(peak_caller, 'seacr'))
    peak_call_func <- peakable::read_seacr
  if (str_detect(peak_caller, 'narrow'))
    peak_call_func <- peakable::read_macs2_narrow
  if (str_detect(peak_caller, 'broad'))
    peak_call_func <- peakable::read_macs2_broad

  message('Peak callers: ', peak_caller)
  message('Get peak ranges ...')

  peak_grl <- bplapply(peak_df$bed_file,
                       peak_call_func, # miss match peak type and func
                       drop_chrM = TRUE,
                       keep_standard_chrom = TRUE,
                       species = species)
  
  names(peak_grl) <- peak_df$sample_id
  
  # remove blacklist: now only works if species is
  if (remove_blacklist) {
    peak_grl <- bplapply(peak_grl,
                         remove_blacklist,
                         blacklist_file = blacklist_file,
                         species = species)
  }
  
  peak_df$number_of_peaks <- sapply(peak_grl, length)

  return(list(df = peak_df, grl = peak_grl))
}


.macs2_bind <- function(macs2_narrow, macs2_broad) {
  # combine macs2_narrow and macs2_broad
  # macs2_narrow and macs2_broad must be the return of peakle_flow
  macs2 <- vector('list', 2)
  macs2$df <- dplyr::bind_rows(macs2_narrow$df, macs2_broad$df)
  macs2$grl <- c(macs2_narrow$grl, macs2_broad$grl)
  return(macs2)
}
