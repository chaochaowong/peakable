#' tools.R
#' tools in this script need to have unit test
#' 


#' @mat matrix (n-by-m)
#' @norm_factor a vector of length m
.normalize_columns <- function(mat, norm_factor) {
  stopifnot(ncol(mat) == length(norm_factor))
  t(t(mat) * norm_factor)
}