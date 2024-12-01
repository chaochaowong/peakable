load_dataset <- function() {
  # return df and grl for a particular peak type

  narrow_files <- list.files(system.file('extdata',
                                         package = 'peakable'),
                             pattern='.narrowPeak$',
                             full.names = TRUE)
  # data frame
  df <- data.frame(cell_line='H1',
                   bed_file = narrow_files,
                   peak_caller = 'mac2_narrow') %>%
    dplyr::mutate(sample_id =
                    stringr::str_replace(basename(bed_file),
                                         '_peaks.narrowPeak', ''),
                  antibody = stringr::str_split(sample_id, '_',
                                                 simplify=TRUE)[, 4],
                  .before=cell_line) %>%
    arrange(antibody)

  # peak ranges
  grl <- lapply(df$bed_file, peakable::read_macs2_narrow,
                drop_chrM = TRUE,
                keep_standard_chrom = TRUE,
                species = 'Homo_sapiens')
  names(grl) <- df$sample_id

  # tidy up
  df <- df %>%
    dplyr::mutate(peak_number = sapply(grl, length)) %>%
    dplyr::select(sample_id, cell_line, antibody, peak_number,
                  peak_caller, bed_file)

  return(list(df=df, grl=grl))
}
