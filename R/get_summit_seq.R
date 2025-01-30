get_summit_seq <- function(gr, bs_genome, extend=200L) {
  # gr is the peaks; only works for macs2 at this moment
  # extract sequence of summit +/- 100
  require(ShortRead)
  summit_ext_200 <- extract_summit_macs2(gr, summit_wid = 200L)
  summit_seq <- getSeq(bs_genome, summit_ext_200)
  names(summit_seq) <- summit_ext_200$name
  summit_seq
}