#' Compute FDR thresholds from genome-wide trans p-values
#'
#' Combines p-values from all chromosomes and uses the `qvalue` package to
#' determine p-value cutoffs at specified FDR levels.
#'
#' @param pval_list Named list of p-value vectors, one per chromosome (from
#'   [aggregate_pvals()]). NAs are removed before analysis.
#' @param fdr_levels Numeric vector of FDR levels to report (default
#'   `c(0.05, 0.1)`).
#'
#' @return A named numeric vector of p-value cutoffs. Names are formatted as
#'   `"fdr0.05"`, `"fdr0.1"`, etc.
#'
#' @export
compute_qvals <- function(pval_list, fdr_levels = c(0.05, 0.1)) {
  allp <- unlist(pval_list)

  # Remove NA, NaN, Inf, and values outside (0, 1]
  n_raw <- length(allp)
  allp  <- allp[is.finite(allp) & allp > 0 & allp <= 1]
  n_dropped <- n_raw - length(allp)
  if (n_dropped > 0) {
    warning(n_dropped, " p-value(s) were non-finite or outside (0, 1] and ",
            "were removed before FDR estimation.")
  }

  if (length(allp) == 0) {
    stop("No valid inter-chromosomal p-values available after filtering. ",
         "Check that aggregate_pvals() produced non-empty results, and that ",
         "at least some genes passed the pearsonR > 0.1 threshold in ",
         "compute_trans_pvals().")
  }

  qobj <- qvalue::qvalue(allp)
  pvalues <- qobj$pvalues
  qvalues <- qobj$qvalues

  cutoffs <- vapply(fdr_levels, function(fdr) {
    flag <- qvalues <= fdr
    if (!any(flag)) return(NA_real_)
    max(pvalues[flag])
  }, numeric(1))

  names(cutoffs) <- paste0("fdr", fdr_levels)
  cutoffs
}
