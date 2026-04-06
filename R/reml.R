# REML and ML variance component estimation functions.
# Ported from func_reml.R (originally based on GCTA-style AI-REML).


# Log of determinant of a matrix (handles negative determinants)
.logdet <- function(p) {
  det_l <- determinant(p, logarithm = TRUE)
  if (det_l$sign[1] == 1) {
    return(as.numeric(det_l$modulus))
  } else {
    return(as.numeric(det_l$modulus) * -1)
  }
}


# Validate inputs for REML/ML functions
.validate_reml_input <- function(A, y, Var, X = NULL) {
  N <- length(y)
  if (length(A) + 1 != length(Var)) {
    stop("Number of variance components does not match number of kinship matrices.")
  }
  for (i in seq_along(A)) {
    if (!isSymmetric(A[[i]])) {
      if (isSymmetric(round(A[[i]], 4))) {
        message("Matrix ", i, " is approximately symmetric, continuing.")
      } else {
        stop("Matrix ", i, " is not symmetric.")
      }
    }
    if (nrow(A[[i]]) != N || ncol(A[[i]]) != N) {
      stop("Matrix ", i, " has wrong dimensions: ", paste(dim(A[[i]]), collapse = "x"))
    }
  }
  if (!is.null(X) && nrow(X) != N) {
    stop("Fixed-effects matrix X has wrong number of rows.")
  }
  TRUE
}


#' Average-Information ML (no fixed effects)
#'
#' Estimates variance components for a linear mixed model using
#' Average-Information REML (GCTA style). Intended for internal use by
#' [run_cvblup()].
#'
#' @param A List of kinship/GRM matrices (N x N each).
#' @param y Phenotype vector (length N).
#' @param Var Initial variance proportions summing to 1 (length = `length(A)+1`).
#' @param verbose Print iteration log (default FALSE).
#' @param CALC_SE Compute standard errors via delta method (default FALSE).
#' @param BOUNDED Stop if any variance component escapes the parameter space.
#'
#' @return A list with `h2` (heritability) and `vc` (variance components).
#'   If `CALC_SE = TRUE`, also includes `se`.
#'
#' @export
aiML <- function(A, y, Var, verbose = FALSE, CALC_SE = FALSE,
                 BOUNDED = FALSE) {

  .validate_reml_input(A, y, Var)

  r <- length(A) + 1
  N <- length(y)
  A[[r]] <- diag(N)

  AI <- matrix(0, ncol = r, nrow = r)
  S  <- matrix(0, ncol = r, nrow = r)
  s  <- matrix(0, ncol = 1, nrow = r)

  l_dif <- 10
  it <- 0
  Var <- var(y) * Var

  # Single EM-REML iteration to initialize
  V <- Reduce("+", mapply(function(a, v) a * v, A, Var, SIMPLIFY = FALSE))
  Vinv <- solve(V)
  P <- Vinv
  logL <- -0.5 * (.logdet(V) + t(y) %*% P %*% y)

  for (i in seq_len(r)) {
    Var[i] <- (Var[i]^2 * t(y) %*% P %*% A[[i]] %*% P %*% y +
                 sum(diag(Var[i] * diag(N) - Var[i]^2 * P %*% A[[i]]))) / N
  }

  V <- Reduce("+", mapply(function(a, v) a * v, A, Var, SIMPLIFY = FALSE))
  Vinv <- solve(V)
  P <- Vinv
  logL <- -0.5 * (.logdet(V) + t(y) %*% P %*% y)

  if (verbose) cat("EM:\t", logL, "\n", sep = "")

  # AI-REML iterations (GCTA style)
  while (it < 100 && (abs(l_dif) >= 1e-4 ||
                      (abs(l_dif) < 1e-2 && l_dif < 0))) {
    it <- it + 1

    for (i in seq_len(r)) {
      for (ii in seq_len(r)) {
        if (i == r && ii == r) {
          AI[r, r] <- t(y) %*% P %*% P %*% P %*% y
        } else if (i == r) {
          AI[r, ii] <- t(y) %*% P %*% P %*% A[[ii]] %*% P %*% y
        } else if (ii == r) {
          AI[i, r] <- t(y) %*% P %*% A[[i]] %*% P %*% P %*% y
        } else {
          AI[i, ii] <- t(y) %*% P %*% A[[i]] %*% P %*% A[[ii]] %*% P %*% y
        }
      }
    }
    AI <- 0.5 * AI

    for (i in seq_len(r)) {
      if (i == r) {
        s[r, 1] <- sum(diag(P)) - (t(y) %*% P %*% P %*% y)
      } else {
        s[i, 1] <- sum(diag(P %*% A[[i]])) -
          (t(y) %*% P %*% A[[i]] %*% P %*% y)
      }
    }
    s <- -0.5 * s

    if (l_dif > 1) {
      Var <- Var + 0.316 * (solve(AI) %*% s)
    } else {
      Var <- Var + solve(AI) %*% s
    }

    V <- Reduce("+", mapply(function(a, v) a * v, A, Var, SIMPLIFY = FALSE))
    Vinv <- solve(V)
    P <- Vinv

    new_logL <- -0.5 * (.logdet(V) + t(y) %*% P %*% y)
    l_dif <- new_logL - logL
    logL <- new_logL

    if (verbose) {
      cat(it, "\t", logL, sep = "")
      for (i in seq_len(r)) cat("\t", Var[i], sep = "")
      cat("\n")
    }

    if (BOUNDED && min(Var / sum(Var)) < 0) {
      if (verbose) message("Variance component escaped parameter space, stopping.")
      break
    }
  }

  if (!CALC_SE) {
    return(list(h2 = Var[1] / sum(Var), vc = Var))
  }

  for (i in seq_len(r)) {
    for (ii in seq_len(r)) {
      S[i, ii] <- sum(diag(P %*% A[[i]] %*% P %*% A[[ii]]))
    }
  }
  S <- 0.5 * S
  Sinv <- solve(S)

  sum_eq <- paste(paste0("x", seq_len(r)), collapse = "+")
  SE.p <- msm::deltamethod(as.formula(paste0("~", sum_eq)), Var, Sinv,
                           ses = TRUE)

  SE.i <- vapply(seq_len(r), function(i) {
    eq <- paste0("~x", i, "/(", sum_eq, ")")
    msm::deltamethod(as.formula(eq), Var, Sinv, ses = TRUE)
  }, numeric(1))

  list(h2 = Var[1] / sum(Var), se = SE.i, vc = Var)
}
