# this script prepare the range and sequence for motif analysis
# range 
#   +/- 50 from summit from macs2 peak calling
#   total 100 bps by extending the max.signal.region (volumn 5) from SEACR

seacr_bed_to_gr <- function(file_name) {
  col_names <- c("chr", "start", "end", "AUC", "max.signal", "max.signal.region", "num")
  gr <- readr::read_delim(file_name, col_names = FALSE) 
  # if the bed file columns are not what SEACR claims to be, then just assume it does follow
  # the "order" of chr, start, end, ...
  # also assume # of column is less thatn 5
  names(gr) <- col_names[1:ncol(gr)]
  gr <- gr %>% tibble::add_column(strand = "*", .name_repair = "minimal")
  gr <- as_granges(gr, seqnames = chr)
}

reframe_summit_region_seacrp <- function(file_name, reframe_width = 200L) {
  # x: must be seacr peak bed files
  require(plyranges)
  require(GenomicRanges)
  peak_gr <- seacr_bed_to_gr(file_name)
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