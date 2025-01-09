remove_blacklist <- function(gr, 
                             blacklist_file = NULL,
                             species = 'Homo_sapiens') {
  # verify species and assign hg38's blacklist if blacklist_file is null
  if (is.null(blacklist_file)) {
    message('Note: peakable only provides blacklist v3.0 for homo_sapiens (hg38)')
    blacklist_file <- .get_blacklist_file(blacklist_file, species) 
  }
  
  # validation: empty or file does not exist
  is_empty <- identical(blacklist_file, '')
  if (is_empty) {
    stop('blacklist_file is empty.')
  }
  
  not_exist <- !file.exists(blacklist_file)
  if (not_exist) {
    stop('blacklist_file: ', blacklist_file, ' does not exists.')
  }
  
  # if blacklist is provided -> check if the format (bed) is valid
  # double check the seqlevel style are the compatible
  # import and make sure it is not empty
  message('Import blacklist file: ', blacklist_file)
  blacklist <- rtracklayer::import.bed(blacklist_file)
  is_compatible <- identical(seqlevelsStyle(blacklist), seqlevelsStyle(gr))
  
  # if not compatible, stop
  if (!is_compatible)
    stop('The sequence levels style of the blacklist and the peak ranges are not compatible.')
  
  gr %>% 
    plyranges::filter_by_non_overlaps(blacklist, minoverlap = 1L)
    
}

.get_blacklist_file <- function(blacklist_file, species) {
  # if null and if species is homo_sapiens
  if (stringr::str_detect(species, 'Homo_sapiens|hg38')) {
    blacklist_file <- 
      system.file('assets', 'blacklists', 'v3.0', 
                  'hg38-blacklist.v3.bed', package='peakable')
    blacklist_file <- file.path(blacklist_file)
  } else {
    blacklist_file <- ''
  }

}

.validate_blacklist_file <- function(blaklist_file) {
  is_empty <- blacklist_file == ''
  not_exist <- !file.exists(blacklist_file)
  
  if (is_empty)
    stop('blacklist_file is empty')
  
  
  if (not_exist)
    stop('blacklist_file: ', blacklist_file, ' does not exists')
  
}

