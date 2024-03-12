.l2_norm <- function(x) {
  sqrt(sum(x^2))
}

cos_similarity <- function(x, y) {
  # sanity check norm
  norms <- .l2_norm(x) * .l2_norm(y)
  if (norms < 1e-25)
    stop('L2 norm of x or y must be greater 0.')
  
  sim <- crossprod(x, y) / norms
  return(as.numeric(sim))
}