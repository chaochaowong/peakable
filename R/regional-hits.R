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

#' @importFrom Rsamtools BamFileList
#' @importFrom GenomicAlignments summarizeOverlaps
# script can be optimized; edit to improve efficiency
peak_read_count <- function(features, sample_df, 
                            spike_in_norm = FALSE) {
  # - sample_df is a data.frame must have sample_id, bam_file, aligned_paired, 
  #   and spike_in_factor (spike_nrom=TRUE) columns
  # - count reads hitting the peakset then normlized by spike_in_norm (if spike_norm=TRUE)
  require(GenomicAlignments)
  # sanity check: file existence; check sample_df columns
  bam_files <- Rsamtools::BamFileList(sample_df$bam_file)
  
  #grl <- GenomicAlignments::GAlignmentPairs(bam_files)
  se <- summarizeOverlaps(features = features,
                          reads = bam_files,
                          ignore.strand = TRUE,
                          singleEnd = FALSE,
                          mode = "Union",
                          fragments=FALSE,
                          inter.feature=FALSE)
  rownames(se) <- paste0('peakname_', 1:length(se))
  colData(se) <- DataFrame(sample_df)
  colnames(se) <- sample_df$sample_id
  
  se$reads <- colSums(assays(se)[['counts']])
  
  if (spike_in_norm) {
    if ('spike_in_factor' %in% names(colData(se))) {
      assays(se)[['spike_in_norm']] <- 
        t(t(assays(se)[['counts']]) * sample_df$spike_in_factor)
      se$reads_spikein_norm <- se$reads * se$spike_in_factor
    }
  }
  
  # if aligned_paried is available
  if ('read_paired' %in% names(colData(se))) {
    se$FRiP <- se$reads / se$read_paired
  }
  
  return(se)

}

#' Construct a hit matrix of peaks against a collection of ranges
#' 
#' @regions: a GRanges object
#' @peaks_grl: a list of GRanges or a GRangesList object representing peaks 
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
  
  regions <- consolidate_peaks(grl)
  
  if (is.null(min_overlap))
    min_overlap <- as.integer(min(width(regions)) / 2)
  
  if (!is.integer(min_overlap))
    min_overlap <- as(min_overLap, 'integer')
  
  mat <- regionhit_per_sample_mat(regions, grl, 
                                  min_overlap=min_overlap)
  return(mat)
}

#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom GenomicRanges countOverlaps
#' @export
regionhit_per_sample_mat <- function(regions, grl, min_overlap = NULL) {
  # same as regional_hits
  # peaks_grl: peaks for each sample
  # regions: a GRanges object, usually merge peaks by a peak caller
  # make sure the seq leveql style of peaks and regions are compatible
  if (is.null(grl))
    names(grl) <- paste0('sample-', seq(1, length(grl)))
  
  regions <- keepStandardChromosomes(regions, pruning.mode="coarse")
  
  if (is.null(min_overlap))
    min_overlap <- as.integer(min(width(regions)) / 2)
  
  
  purrr::map_dfc(grl, function(pks) {
    pks <- keepStandardChromosomes(pks, pruning.mode="coarse")
    #seqlevelsStyle(pks) <- 'Ensembl'
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

#' @param mat, sample_info, n_pcs
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join
.getPCA <- function(mat, sample_info=NULL, n_pcs=2) {
  pca <- prcomp(mat, scale. = TRUE)
  pcs <- as.data.frame(pca$rotation[, 1:n_pcs]) 
  
  if (!is.null(sample_info) & 'sample_id' %in% names(sample_info)) {
    pcs <- pcs %>%
      tibble::rownames_to_column(var='sample_id') %>%
      dplyr::left_join(sample_info, by='sample_id')
  }
  return(pcs)
}

