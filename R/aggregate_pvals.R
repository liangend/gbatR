#' Aggregate inter-chromosomal trans p-values for one chromosome
#'
#' Combines p-values from all subsets of a chromosome, removes genes with
#' invalid heritability estimates, and filters out same-chromosome
#' (potentially cis) associations.
#'
#' @param cor_results_list Named list of p-value matrices from
#'   [compute_trans_pvals()].
#' @param h2g_results_list Named list of h2g data frames from [run_cvblup()].
#'   Must align with `cor_results_list`.
#' @param gene_pos Data frame with columns for gene ID (specified by
#'   `gene_col`) and `chr`. Used to identify which genes are on the same
#'   chromosome as each predictor and should be excluded from trans p-values.
#' @param chr Chromosome number of the predictor genes (integer or character).
#' @param gene_col Name of the gene ID column in `gene_pos` (default `"gene"`).
#'
#' @return A numeric vector of inter-chromosomal trans p-values.
#'
#' @export
aggregate_pvals <- function(cor_results_list, h2g_results_list,
                            gene_pos, chr, gene_col = "gene") {
  chrom_str <- paste0("chr", chr)
  allp      <- NULL

  for (i in seq_along(cor_results_list)) {
    dat <- cor_results_list[[i]]
    h2  <- h2g_results_list[[i]]

    if (is.null(dat) || ncol(dat) == 0) next

    # Remove genes with invalid h2g (h2g column always named "h2g")
    # The gene ID column name in h2 matches gene_col
    neg_flag  <- which(h2$h2g < 0 | h2$h2g > 1)
    neg_genes <- as.character(h2[[gene_col]][neg_flag])
    gene_names <- colnames(dat)
    gene_flag  <- which(gene_names %in% neg_genes)
    if (length(gene_flag) > 0) {
      dat        <- dat[, -gene_flag, drop = FALSE]
      gene_names <- colnames(dat)
    }

    # For each predictor gene: drop rows on the same chromosome
    gp_ids <- gene_pos[[gene_col]]
    for (nc in seq_len(ncol(dat))) {
      gene <- gene_names[nc]
      flag <- which(gp_ids == gene)
      if (length(flag) == 0) next
      flag_chr <- gene_pos$chr[flag[1]]

      same_chr_genes <- gp_ids[gene_pos$chr == flag_chr]
      rm_flag <- which(rownames(dat) %in% same_chr_genes)
      row_idx <- setdiff(seq_len(nrow(dat)), rm_flag)
      allp <- c(allp, as.numeric(dat[row_idx, nc]))
    }
  }

  allp
}
