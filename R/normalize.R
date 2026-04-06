# Normalize a raw count matrix into expression values
#
# Starting from a raw read count matrix, removes sex/mitochondrial chromosome
# genes, filters lowly expressed genes, computes RPKM, and applies two rounds
# of quantile normalization.
#
# @param count_matrix Numeric matrix of raw counts, genes x samples.
#   Row names must be gene IDs matching `gene_info$gene`.
# @param gene_info Data frame with columns: `gene`, `chr`, `start`, `end`,
#   `length`. Chromosome values should be in "chr1" format.
# @param min_cpm Minimum CPM threshold for expression filter (default 0.5).
# @param min_sample_frac Minimum fraction of samples that must meet `min_cpm`
#   for a gene to be retained (default 0.5).
#
# @return A list with:
#   {expr}{Normalized expression matrix, samples x genes.}
#   {gene_info}{Filtered gene annotation data frame, matching columns
#     of `expr`.}

normalize_counts <- function(count_matrix, gene_info,
                             min_cpm = 0.5, min_sample_frac = 0.5) {

  # Align count matrix rows with gene_info
  common_genes <- intersect(rownames(count_matrix), gene_info$gene)
  if (length(common_genes) == 0) {
    stop("No genes match between count_matrix rownames and gene_info$gene.")
  }
  count_matrix <- count_matrix[common_genes, , drop = FALSE]
  gene_info <- gene_info[match(common_genes, gene_info$gene), ]

  # Remove sex and mitochondrial chromosomes
  rm_flag <- gene_info$chr %in% c("chrX", "chrY", "chrM")
  if (any(rm_flag)) {
    count_matrix <- count_matrix[!rm_flag, , drop = FALSE]
    gene_info <- gene_info[!rm_flag, ]
  }

  N <- ncol(count_matrix)

  # CPM using total counts per sample as the library size
  total_count <- colSums(count_matrix)
  cpm <- sweep(count_matrix, 2, total_count / 1e6, "/")

  # Keep genes expressed (CPM >= min_cpm) in at least min_sample_frac samples
  ex_count <- rowSums(cpm >= min_cpm)
  keep_flag <- ex_count >= N * min_sample_frac
  count_matrix <- count_matrix[keep_flag, , drop = FALSE]
  gene_info <- gene_info[keep_flag, ]

  # RPKM: normalize by gene length, then by total RPK per sample
  gene_length <- abs(gene_info$length)
  rpk <- sweep(count_matrix, 1, gene_length, "/")
  total_rpk <- colSums(rpk)
  rpkm <- sweep(rpk, 2, total_rpk / 1e6, "/")

  # Round 1: quantile normalize across samples (genes x samples -> equal
  # distributions per sample)
  qn <- limma::normalizeQuantiles(as.matrix(rpkm))

  # Round 2: rank normalize across genes (samples x genes -> each gene
  # mapped to standard normal)
  expr <- apply(qn, 1, RNOmni::RankNorm)
  # expr is now samples x genes

  rownames(expr) <- colnames(count_matrix)
  colnames(expr) <- gene_info$gene

  list(expr = expr, gene_info = gene_info)
}
