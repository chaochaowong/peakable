#' import MACS2 narrowPeak and broadPeak bed file
#' 

read_macs2_narrow <- function(file, ...) {
  # sanity check
  stopifnot(length(file) == 1)
  stopifnot(file.exists(file)) 
  
  # define extra columns
  extraCols_narrowPeak <- c(signalValue = "numeric", 
                            pValue = "numeric",
                            qValue = "numeric", 
                            peak = "integer")
  
  rtracklayer::import(file, format = 'BED', 
                      extraCols = extraCols_narrowPeak, ...) 
}

read_macs2_broad <- function(file, ...) {
  # sanity check
  stopifnot(length(file) == 1)
  stopifnot(file.exists(file)) 
  
  # define extra columns
  extraCols_broadPeak <- c(signalValue = "numeric", 
                           pValue = "numeric",
                           qValue = "numeric")
  
  rtracklayer::import(file, format = 'BED',
                      extraCols = extraCols_broadPeak, ...) 
}
