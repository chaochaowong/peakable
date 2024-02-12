#' @usage consolidate_peaks(grl)
#' @param grl
#' @grl A list of GRanges or a GRangeList object
#' @return GRanges
#'
#' @importFeom GenomicRanges reduce
#' @export
consolidate_peaks <- function(grl) {
  # should update to as S3 and S4 method
  # consolidate peaks from a list of GRanges or GRangesList
  # keep standard chromosomes
  # keep clean: remove mcols and only standard chromosome
  gr <- unlist(as(grl, "GRangesList"))
  gr <- keepStandardChromosomes(gr, pruning.mode='coarse')
  mcols(gr) <- NULL
  GenomicRanges::reduce(gr, ignore.strand=TRUE)
}


#'
#' consolidated_peak_hits
#' Construct a hit matrix of peaks against a collection of ranges
#' 
#' @regions: a GRanges object
#' @grl: a list of GRanges or a GRangesList object representing peaks 
#' @min_overlap: minimal overal between peaks and regions
#' 
#' @return a matrix
#' @importFram purrr map_df
#' @example
#' @export
consolidated_peak_hits <- function(grl, min_overlap = NULL) {
  
  # grl: peaks for each sample
  # make sure the seq level style of peaks and regions are compatible
  
  # need unit test:
  # the col names of mat matches the names pf peaks_grl
  # the seq levels of peaks_grl must match each other
  
  # check names of grl
  if (is.null(grl))
    stop('The names of grl is NULL. The names should be the sample ID and must be given.')
  
  #' consolidate all the peaks from grl; only keep standard chromosomes
  regions <- consolidate_peaks(grl)
  
  if (is.null(min_overlap))
    min_overlap <- as.integer(min(width(regions)) / 2)
  
  if (!is.integer(min_overlap))
    min_overlap <- as(min_overLap, 'integer')
  
  mat <- regionhit_per_sample_mat(regions, grl, 
                                  min_overlap=min_overlap)
  return(mat)
}


#' regionhit_per_sample_mat
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom GenomicRanges countOverlaps
#' @importFrom purrr map_dfc
#' @export
regionhit_per_sample_mat <- function(regions, grl, min_overlap = NULL) {
  # same as regional_hits
  # peaks_grl: peaks for each sample
  # regions: a GRanges object, usually merge peaks by a peak caller
  # make sure the seq leveql style of peaks and regions are compatible
  if (is.null(names(grl)))
    names(grl) <- paste0('sample-', seq(1, length(grl)))
  
  regions <- keepStandardChromosomes(regions, pruning.mode="coarse")
  
  if (is.null(min_overlap))
    min_overlap <- as.integer(min(width(regions)) / 2)
  
  
  purrr::map_dfc(grl, function(pks) {
    pks <- keepStandardChromosomes(pks, pruning.mode="coarse")
    cnt <- countOverlaps(query=regions, subject=pks, 
                         ignore.strand=TRUE,
                         type = 'any',
                         minoverlap = min_overlap)
    cnt <- if_else(cnt > 0, 1, 0) # convert number greater than 0 to 1
  })
}

promoterhit_per_sample_mat <- function(ensdb, peaks_grl,
                                       upstream=3000,
                                       downstream=300) {
  #' for Ensembl only at this moment
  #' construct promoter-hit-per-sample matrix (0,1)
  #' Input: 
  #'     ensdb: EnsDb or TxDb
  require(purrr)
  require(GenomiceRanges)

  around_tss <- promoters(ensdb, 
                          upstream = upsream,
                          downstream = downstream, 
                          use.names=TRUE, 
                          columns = c("gene_id", "gene_name", "gene_biotype"))
  around_tss <- keepStandardChromosomes(around_tss, pruning.mode="coarse")
  
  tmp <- map_df(peaks_grl, function(pks) {
    pks <- keepStandardChromosomes(pks, pruning.mode="coarse")
    seqlevelsStyle(pks) <- 'Ensembl'
    cnt <- countOverlaps(query=around_tss, subject=pks, ignore.strand=TRUE,
                         minoverlap = 100L)
    if_else(cnt > 0, 1, 0) # convert number greater than 0 to 1
  }) %>% 
    add_column(gene_id = around_tss$gene_id, gene_name= around_tss$gene_name)
  tx_per_sample <- tmp[rowSums(tmp[, 1:(ncol(tmp)-2)]) > 0, ] 
}

genehit_per_sample_mat <- function(ensdb, peaks_grl) {
  gene_rng <- genes(ensdb, columns = c(listColumns(ensdb, 'gene'), 'gene_name'))
  # only keep the standard chromosomes
  gene_rng <- keepStandardChromosomes(gene_rng, pruning.mode="coarse")
  # extend 3000 bps from 5p: anchor at 3p and extend 3000 bps
  gene_rng <- plyranges::stretch(anchor_3p(gene_rng), 3000)
  
  tmp <- map_df(peaks_grl, function(pks) {
    pks <- keepStandardChromosomes(pks, pruning.mode="coarse")
    seqlevelsStyle(pks) <- 'Ensembl'
    cnt <- countOverlaps(query=gene_rng, subject=pks, ignore.strand=TRUE,
                         minoverlap = 100L)
    if_else(cnt > 0, 1, 0)
  }) %>% 
    add_column(gene_id = gene_rng$gene_id, gene_name= gene_rng$gene_name)
  gene_per_sample = tmp[rowSums(tmp[, 1:(ncol(tmp)-2)]) > 0, ] 
}



