# script can be optimized; edit to improve efficiency
peakCoverageMatrixRSE <- function(peakset_gr, sample_df, spike_in_norm = TRUE,
                               spike_in_factor = NULL) {
  # - sample_df is a data.frame must have sample_id, bam_file, aligned_paired, 
  #   and spike_in_factor (spike_nrom=TRUE) columns
  # - count reads hitting the peakset then normlized by spike_in_norm (if spike_norm=TRUE)
  require(GenomicAlignments)
  # sanity check: file existence; check sample_df columns
  bam_files <- BamFileList(sample_df$bam_file)
  
  #grl <- GenomicAlignments::GAlignmentPairs(bam_files)
  se <- summarizeOverlaps(features = peakset_gr,
                          reads = bam_files,
                          ignore.strand = TRUE,
                          singleEnd = FALSE,
                          mode = "Union",
                          fragments=FALSE,
                          inter.feature=FALSE)
  assays(se)[['spikein_norm']] <- t(t(assays(se)[['counts']]) * sample_df$spikein_factor)
  
  rownames(se) <- paste0('peakname_', 1:length(se))
  colData(se) <- DataFrame(sample_df)
  se$reads <- colSums(assays(se)[['counts']])
  se$reads_spikein_norm <- se$reads * se$spikein_factor
  # if aligned_paried is available
  se$FRiP <- se$reads / se$aligned_paired
  
  return(se)

}

regionhit_per_sample_mat <- function(regions, peaks_grl) {
  # peaks_grl: peaks for each sample
  # regions: usually merge peaks by a peak caller
  require(purrr)
  require(GenomicRanges)
  regions <- keepStandardChromosomes(regions, pruning.mode="coarse")
  min_overlap <- as.integer(min(width(regions)) / 2)
  map_df(peaks_grl, function(pks) {
    pks <- keepStandardChromosomes(pks, pruning.mode="coarse")
    #seqlevelsStyle(pks) <- 'Ensembl'
    cnt <- countOverlaps(query=regions, subject=pks, ignore.strand=TRUE,
                         minoverlap = min_overlap)
    cnt <- if_else(cnt > 0, 1, 0) # convert number greater than 0 to 1
    return(cnt)
  })
}

promoterhit_per_sample_mat <- function(ensdb, peaks_grl) {
  #' construct promoter-hit-per-sample matrix (0,1)
  #' Input: 
  #'     ensdb: EnsDb or TxDb
  require(purrr)
  require(GenomiceRanges)

  around_tss <- promoters(ensdb, upstream = 3000, downstream=3000, 
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

.pca_pc1_pc2 <- function(mat, tb_to_join, n_pcs=2) {
  pca <- prcomp(mat, scale. = TRUE)
  as.data.frame(pca$rotation[, 1:n_pcs]) %>%
    rownames_to_column(var='callpeaks') %>%
    left_join(tb_to_join, by='callpeaks')
}

# for each gene or promoter (depending on the mat), value is 1 if both duplicates hitting on it
# region-per-sample: combine duplicates—if both 2, then TRUE, otherwise FALSE
.comb_duplicates_matrix <- function(hit_mat, duplicates) {
  duplicate <- c("PDX_2_sample_9", "PDX_2_sample_10",
                "RCH", "REH")
  tmp <- map_dfc(duplicate, function(x) {
    idx <- grepl(x, names(hit_mat))
    rowSums(hit_mat[, idx]) == 2
    #hit_mat[, idx] %>% mutate_all(as.logical) %>% 
     # rowSums() 
  }) %>%
    setNames(duplicate) %>%
    add_column(gene_id = hit_mat$gene_id,
               gene_name = hit_mat$gene_name)
}