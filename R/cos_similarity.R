#' cos_similarity
#'
#' @param x GRanges
#' @param y GRanges
#' @export
cos_similarity <- function(x, y) {
  # Get corpus peaks and the hit matrix
  hits <- consolidated_peak_hits(list(x, y))
  x <- hits[, 1, drop=TRUE]
  y <- hits[, 2, drop=TRUE]

  # L2 norm
  norms <- .L2_norm(x) * .L2_norm(y)

  # sanity check
  if (norms < 1e-25)
    stop('L2 norm of x or y must be greater 0.')

  # cosine(theta)
  sim <- crossprod(x, y) / norms

  return(as.numeric(sim))
}

.L2_norm <- function(x) {
  sqrt(sum(x^2))
}
