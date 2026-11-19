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
#' @param num_sv Number of supervised SVs to estimate and use (default 20).
#' @param covariates Optional numeric sample-by-covariate matrix. Row names
#'   are sample IDs and are matched to expression before testing.
#'
#' @return A matrix of trans p-values. Rows = all expression genes, columns =
#'   tested predictor genes (those with sufficient cis signal). Row names are
#'   expression gene IDs; column names are predictor gene IDs.
#'
#' @export
compute_trans_pvals <- function(expr, blup_matrix, pearson_r, num_sv = 20,
                                covariates = NULL) {
  if (length(num_sv) != 1L || !is.finite(num_sv) ||
      num_sv != floor(num_sv) || num_sv < 0) {
    stop("num_sv must be a non-negative integer.")
  }
  covariates <- .align_covariates(covariates, expr)
  if (!is.null(rownames(blup_matrix)) || !is.null(covariates)) {
    blup_matrix <- .align_samples(blup_matrix, expr, "BLUPs")
  } else if (nrow(blup_matrix) != nrow(expr)) {
    stop("BLUP and expression matrices must contain the same samples.")
  }
  design <- cbind(intercept = rep(1, nrow(expr)), covariates)
  if (num_sv > 0L &&
      (ncol(expr) < num_sv || nrow(expr) <= ncol(design) + num_sv + 1L)) {
    stop("Too few expression genes or samples to estimate ", num_sv,
         " supervised SVs with the supplied covariates.")
  }
  if (qr(design)$rank < ncol(design)) {
    stop("Covariates must be linearly independent and exclude an intercept.")
  }
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
      g_sva <- if (num_sv > 0L) {
        .make_smartsva(expr, num_sv, g_blup, covariates)
      } else {
        matrix(numeric(0), nrow = nrow(expr), ncol = 0L)
      }

      # Standardize BLUP
      g_norm <- (g_blup - mean(g_blup)) / blup_sd

      if (is.null(dim(g_sva)) || ncol(g_sva) < num_sv ||
          any(!is.finite(g_sva))) {
        stop("SmartSVA did not return a finite surrogate-variable matrix.")
      }
      adjustment <- cbind(covariates, g_sva[, seq_len(num_sv), drop = FALSE])
      if (ncol(adjustment) > 0L) {
        rlm <- stats::lm(expr ~ adjustment)
        rQN <- phenix::quantnorm(as.matrix(stats::resid(rlm)))
      } else {
        rQN <- expr
      }

      # Keep known covariates in the final test after quantile normalization.
      test <- if (is.null(covariates)) {
        stats::lm(rQN ~ g_norm)
      } else {
        stats::lm(rQN ~ covariates + g_norm)
      }
      pvals_col <- sapply(summary(test),
                          function(x) if ("g_norm" %in% rownames(x$coefficients)) {
                            x$coefficients["g_norm", 4]
                          } else {
                            NA_real_
                          })

      all_cor <- cbind(all_cor, pvals_col)
      coln <- c(coln, gene)
    }, error = function(e) {
      warning("Trans test failed for ", gene, ": ", conditionMessage(e),
              call. = FALSE)
    })
  }

  if (is.null(all_cor)) {
    return(matrix(nrow = ncol(expr), ncol = 0,
                  dimnames = list(colnames(expr), NULL)))
  }

  rownames(all_cor) <- colnames(expr)
  colnames(all_cor) <- coln
  all_cor
}
