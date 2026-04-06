#' Extract significant trans-eQTL regulator-target pairs
#'
#' Filters trans p-value matrices at a given significance cutoff and returns
#' one row per significant regulator-target gene pair. In the p-value matrix,
#' columns are regulator genes (those with BLUP predictions from cis-SNPs) and
#' rows are target genes (all expression genes). Pairs on the same chromosome
#' within `cis_exclusion_window` bp are excluded.
#'
#' @param cor_results_list Named list of trans p-value matrices from
#'   [compute_trans_pvals()]. Each matrix: rows = target genes, columns =
#'   regulator genes.
#' @param gene_pos Data frame with columns for gene ID (specified by
#'   `gene_col`), `chr`, `start`, `end`. Must contain entries for both
#'   regulator and target genes.
#' @param cutoff P-value cutoff (e.g. from [compute_qvals()]).
#' @param gene_col Name of the gene ID column in `gene_pos` (default `"gene"`).
#' @param cis_exclusion_window Window (bp) used to exclude regulator-target
#'   pairs on the same chromosome that may reflect cis effects
#'   (default 10,000,000 bp = 10 Mb).
#' @param pseudogenes Optional character vector of gene IDs to exclude from
#'   both regulators and targets. Pass `NULL` to skip.
#'
#' @return A data frame with one row per significant regulator-target pair,
#'   ordered by regulator chromosome and start position, with columns:
#'   \item{regulator_gene}{Regulator gene ID.}
#'   \item{regulator_chr}{Regulator gene chromosome.}
#'   \item{regulator_start}{Regulator gene start position.}
#'   \item{regulator_end}{Regulator gene end position.}
#'   \item{target_gene}{Target gene ID.}
#'   \item{target_chr}{Target gene chromosome.}
#'   \item{target_start}{Target gene start position.}
#'   \item{target_end}{Target gene end position.}
#'   \item{pval}{GBAT trans p-value for this pair.}
#'
#' @export
get_sig_genes <- function(cor_results_list, gene_pos,
                          cutoff,
                          gene_col = "gene",
                          cis_exclusion_window = 1e7,
                          pseudogenes = NULL) {

  results <- NULL

  for (sub_name in names(cor_results_list)) {
    dat <- cor_results_list[[sub_name]]

    if (is.null(dat) || ncol(dat) == 0) next

    dat <- as.matrix(dat)

    # Iterate over columns (regulator genes)
    for (j in seq_len(ncol(dat))) {
      reg_gene <- colnames(dat)[j]

      if (!is.null(pseudogenes) && reg_gene %in% pseudogenes) next

      reg_info <- gene_pos[gene_pos[[gene_col]] == reg_gene, ][1L, ]
      if (is.na(reg_info[[gene_col]])) next

      # Find significant rows (target genes) for this regulator
      sig_rows <- which(dat[, j] < cutoff)
      if (length(sig_rows) == 0) next

      for (i in sig_rows) {
        target_gene <- rownames(dat)[i]

        if (!is.null(pseudogenes) && target_gene %in% pseudogenes) next

        target_info <- gene_pos[gene_pos[[gene_col]] == target_gene, ][1L, ]
        if (is.na(target_info[[gene_col]])) next

        # Exclude pairs on the same chromosome within cis_exclusion_window
        # Use end position if available, otherwise fall back to start
        if (reg_info$chr == target_info$chr) {
          reg_end    <- if (!is.null(reg_info$end))    reg_info$end    else reg_info$start
          target_end <- if (!is.null(target_info$end)) target_info$end else target_info$start
          overlap <- (reg_info$start > target_info$start - cis_exclusion_window &
                        reg_info$start < target_end         + cis_exclusion_window) |
                     (reg_end         > target_info$start - cis_exclusion_window &
                        reg_end         < target_end         + cis_exclusion_window)
          if (isTRUE(any(overlap))) next
        }

        results <- rbind(results, data.frame(
          regulator_gene  = reg_gene,
          regulator_chr   = reg_info$chr,
          regulator_start = reg_info$start,
          regulator_end   = if (!is.null(reg_info$end))    reg_info$end    else NA_integer_,
          target_gene     = target_gene,
          target_chr      = target_info$chr,
          target_start    = target_info$start,
          target_end      = if (!is.null(target_info$end)) target_info$end else NA_integer_,
          pval            = dat[i, j],
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # Order by regulator chromosome (numeric) then regulator start position
  if (!is.null(results) && nrow(results) > 0) {
    ord <- order(
      as.integer(sub("chr", "", results$regulator_chr)),
      results$regulator_start
    )
    results <- results[ord, ]
    rownames(results) <- NULL
  }

  results
}
