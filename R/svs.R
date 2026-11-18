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


# Match a sample-indexed matrix to expression, allowing additional input rows.
.align_samples <- function(x, expr, label) {
  ids <- rownames(expr)
  input_ids <- rownames(x)
  valid_ids <- function(z) {
    !is.null(z) && !anyNA(z) && all(nzchar(z)) && !anyDuplicated(z)
  }
  if (!valid_ids(ids) || !valid_ids(input_ids)) {
    stop("Expression and ", label, " must have unique, non-empty sample IDs as row names.")
  }
  index <- match(ids, input_ids)
  if (anyNA(index)) {
    stop("Missing sample IDs in ", label, ": ",
         paste(head(ids[is.na(index)]), collapse = ", "))
  }
  x[index, , drop = FALSE]
}

.align_covariates <- function(covariates, expr) {
  if (is.null(covariates)) return(NULL)
  covariates <- .align_samples(as.matrix(covariates), expr, "covariates")
  if (!is.numeric(covariates) || any(!is.finite(covariates))) {
    stop("Covariates must be numeric and finite for all expression samples.")
  }
  covariates
}

# Estimate 100 supervised SVs with known covariates in both designs.
.make_smartsva <- function(expr, K, X, covariates = NULL) {
  mod0 <- cbind(intercept = rep(1, nrow(expr)), covariates)
  mod <- cbind(mod0, predictor = X)
  SmartSVA::smartsva.cpp(dat = t(expr), mod = mod, n.sv = K, mod0 = mod0)$sv
}
