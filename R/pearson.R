#' Compute Pearson correlation between BLUP predictions and expression
#'
#' For each gene's BLUP prediction, computes the Pearson correlation against
#' the corresponding column in the expression matrix. This is used as a
#' quality-control metric: a BLUP with very low correlation to its own
#' gene's expression indicates weak cis-genetic signal.
#'
#' @param blup_matrix N x G matrix of BLUP predictions from [run_cvblup()].
#'   Column names are gene IDs.
#' @param expr Normalized expression matrix, samples x genes. Column names
#'   must include the genes in `blup_matrix`.
#'
#' @return A data frame with columns:
#'   \item{gene}{Gene ID.}
#'   \item{pearsonR}{Pearson correlation between the gene's BLUP and its
#'     observed expression.}
#'
#' @export
compute_pearson_r <- function(blup_matrix, expr) {
  gene_names <- colnames(expr)
  pearson_r <- vapply(seq_len(ncol(blup_matrix)), function(i) {
    gene <- colnames(blup_matrix)[i]
    ex_flag <- which(gene_names == gene)
    if (length(ex_flag) == 0) return(NA_real_)
    stats::cor(blup_matrix[, i], expr[, ex_flag])
  }, numeric(1))

  data.frame(
    gene     = colnames(blup_matrix),
    pearsonR = pearson_r,
    stringsAsFactors = FALSE
  )
}
