read_fire_peaks <- function(file, col_names = NULL,
                            FDR_threshold = 0.01,
                            genome_info = NULL) {
  #args <- norm_args_reader(genome_info)
    
  df <- read_tsv(file) %>%
    dplyr::rename(seqnames = `#chrom`) %>%
    dplyr::mutate(strand = '*') %>%
    dplyr::filter(FDR < FDR_threshold)
  
  as_granges(df)
}