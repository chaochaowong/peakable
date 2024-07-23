# this script prepare the range and sequence for motif analysis
# range
#   +/- 50 from summit from macs2 peak calling
#   total 100 bps by extending the max.signal.region (volumn 5) from SEACR

reframe_summit_region_seacrp <- function(file_name, reframe_width = 200L) {
  # x: must be seacr peak bed files
  require(plyranges)
  require(GenomicRanges)
  peak_gr <- read_seacr(file_name)
  summit_gr <- GRanges(peak_gr$max.signal.region)
  summit_gr <- plyranges::mutate(anchor_center(summit_gr),
                                 width = reframe_width) %>%
    mutate(score = mcols(peak_gr)$max.signal) %>%
    arrange(desc(score)) %>%
    mutate(name = paste0('peak_', 1:length(.)))
}

reframe_summit_to_bed <- function(file_name, reframe_width = 200L,
                                  path) {
  require(rtracklayer)
  require(plyranges)
  require(GenomicRanges)
  summit_gr <- reframe_summit_region_seacrp(file_name, reframe_width = 200L)
  rtracklayer::export(summit_gr, path, format='narrowPeak')
  summit_gr
}

reframe_summit_to_fa <- function(file_name, reframe_width = 200L,
                                 top_n = 1000L, fa_path) {
  require(rtracklayer)
  require(plyranges)
  require(GenomicRanges)
  require(BSgenome.Hsapiens.UCSC.hg38)
  require(Biostrings)
  genome <- BSgenome.Hsapiens.UCSC.hg38

  summit_gr <- reframe_summit_region_seacrp(file_name,
                                            reframe_width = reframe_width)

  if (length(summit_gr) > top_n) {
    summit_gr <- summit_gr[1:top_n]
  }

  summit_seq <- getSeq(genome, summit_gr)
  names(summit_seq) <- summit_gr$name

  writeXStringSet(summit_seq, file = fa_path, format = "fasta")
  invisible()
}
