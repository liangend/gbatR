#' Compute surrogate variables from expression data
#'
#' Runs SVA (Surrogate Variable Analysis) using the two-step method to
#' estimate latent factors that capture unwanted variation in expression data.
#'
#' @param expr Normalized expression matrix, samples x genes.
#' @param num_sv Number of surrogate variables to estimate (default 20).
#'
#' @return A samples x `num_sv` matrix of surrogate variables, or `NA` if
#'   all estimated SVs are zero.
#'
#' @export
compute_svs <- function(expr, num_sv = 20) {
  mod0 <- model.matrix(~1, data = data.frame(expr))
  out <- sva::sva(dat = t(expr), mod = mod0, n.sv = num_sv,
                  method = "two-step")
  svs <- out$sv
  if (all(svs == 0) && is.null(dim(svs))) {
    return(NA)
  }
  svs
}


# Internal: run SmartSVA with an optional continuous covariate X.
# Used inside compute_trans_pvals() to re-estimate SVs per BLUP column.
.make_smartsva <- function(expr, K, X) {
  mod0 <- model.matrix(~1, data = data.frame(expr))
  if (!missing(X)) {
    mod <- model.matrix(~1 + X, data = data.frame(expr))
    out <- SmartSVA::smartsva.cpp(dat = t(expr), mod = mod,
                                  n.sv = K, mod0 = mod0)
  } else {
    out <- SmartSVA::smartsva.cpp(dat = t(expr), mod = mod0, n.sv = K)
  }
  svs <- out$sv
  if (all(svs == 0) && is.null(dim(svs))) {
    return(NA)
  }
  svs
}
