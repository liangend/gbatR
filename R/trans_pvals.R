#' Compute trans-association p-values via SmartSVA residualization
#'
#' For each gene with a BLUP prediction (and sufficient cis-heritability),
#' re-estimates surrogate variables using SmartSVA with the BLUP as a
#' covariate, residualizes the expression matrix, then regresses residualized
#' expression of all genes on the standardized BLUP to obtain trans p-values.
#'
#' @param expr Normalized expression matrix, samples x genes.
#' @param blup_matrix N x G matrix of BLUP predictions from [run_cvblup()].
#' @param pearson_r Data frame from [compute_pearson_r()] with columns
#'   `gene` and `pearsonR`. Only genes with `pearsonR > 0.1` are tested.
#' @param num_sv Number of surrogate variables for SmartSVA (default 20).
#'
#' @return A matrix of trans p-values. Rows = all expression genes, columns =
#'   tested predictor genes (those with sufficient cis signal). Row names are
#'   expression gene IDs; column names are predictor gene IDs.
#'
#' @export
compute_trans_pvals <- function(expr, blup_matrix, pearson_r, num_sv = 20) {
  all_cor <- NULL
  coln <- NULL

  for (i in seq_len(ncol(blup_matrix))) {
    gene <- colnames(blup_matrix)[i]

    # Skip if BLUP is all NA or Pearson R is too low
    if (all(is.na(blup_matrix[, i]))) next
    pr_i <- pearson_r$pearsonR[pearson_r$gene == gene]
    if (length(pr_i) == 0 || is.na(pr_i) || pr_i <= 0.1) next

    tryCatch({
      g_blup <- blup_matrix[, i]

      # Skip if BLUP has zero or near-zero variance (cannot standardize)
      blup_sd <- stats::sd(g_blup)
      if (is.na(blup_sd) || blup_sd < 1e-10) next

      # Re-estimate SVs using BLUP as a covariate (SmartSVA)
      g_sva <- .make_smartsva(expr, num_sv, g_blup)

      # Standardize BLUP
      g_norm <- (g_blup - mean(g_blup)) / blup_sd

      if (!any(is.na(g_sva))) {
        # Residualize expression against SVA (first num_sv columns)
        n_sv_use <- min(num_sv, ncol(g_sva))
        rlm <- stats::lm(expr ~ g_sva[, seq_len(n_sv_use)])
        rQN <- phenix::quantnorm(as.matrix(stats::resid(rlm)))
      } else {
        rQN <- expr
      }

      # Regress residualized expression on standardized BLUP
      test <- stats::lm(rQN ~ g_norm)
      pvals_col <- sapply(summary(test),
                          function(x) x$coefficients[2, 4])

      all_cor <- cbind(all_cor, pvals_col)
      coln <- c(coln, gene)
    }, error = function(e) {})
  }

  if (is.null(all_cor)) {
    return(matrix(nrow = nrow(expr), ncol = 0))
  }

  rownames(all_cor) <- colnames(expr)
  colnames(all_cor) <- coln
  all_cor
}
