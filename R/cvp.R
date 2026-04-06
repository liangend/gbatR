# Cross-validated BLUP functions.
# Ported from cvp.r in the original GBAT pipeline.
#
# Key function: .getBlups()
#
# Computes leave-one-out cross-validated BLUP predictions for a linear
# mixed model:  y = X*beta + Z*b + epsilon
#
# Fixed effects (beta) are also cross-validated via LOO, then given those
# residuals, the genetic BLUPs (gBlupLoo) are estimated by LOO.
# sig2g and sig2e are estimated once using all data.


# Internal helper: n-by-p matrix of ones
.ones <- function(n, p) matrix(1, n, p)


# Internal: LOO bBLUPs (SNP-level effect estimates, m-by-n matrix).
# Used only when return_bBlupsLoo = TRUE inside .getBlups().
.getBBlupsLoo <- function(y, Z, Vi, yfeLoo, gBlupLoo, sig2g) {
  y      <- as.numeric(y)
  Z      <- as.matrix(Z)
  n      <- nrow(Z)
  m      <- ncol(Z)
  yTilde <- as.numeric(y - yfeLoo)
  Y      <- as.matrix(diag(yTilde)) %*% matrix(1, n, n)
  diag(Y) <- gBlupLoo
  sig2g / m * t(Z) %*% Vi %*% Y
}


# Compute cross-validated BLUPs (LOO scheme).
#
# y   : numeric phenotype vector, length n
# X   : n-by-p covariate matrix (intercept added internally), or NA
# Z   : n-by-m matrix of scaled, centred genotypes
# sig2g, sig2e : variance components; supply h2 instead if unknown
# h2  : heritability; used to derive sig2g/sig2e when those are NA
# return_bBlupsLoo : if TRUE also return m-by-n SNP BLUP LOO matrix
# return_intermediates : if TRUE also return V, Vi, H, sig2g, sig2e
#
# Returns a list including gBlupLoo (the LOO genetic predictions, length n).
.getBlups <- function(y, X = NA, Z, sig2g = NA, sig2e = NA, h2 = NA,
                      return_bBlupsLoo = FALSE,
                      return_intermediates = FALSE) {

  y <- as.numeric(y)
  Z <- as.matrix(Z)
  n <- nrow(Z)
  m <- ncol(Z)
  K <- 1 / m * tcrossprod(Z)

  # Build design matrix with intercept
  if (all(is.na(X))) {
    X <- as.matrix(data.frame(rep(1, n)))
  } else {
    X <- as.matrix(data.frame(cbind(rep(1, n), X)))
  }

  # Estimate sig2g / sig2e from h2 if not supplied
  if (is.na(sig2g) || is.na(sig2e)) {
    betaOLS  <- solve(crossprod(X)) %*% (t(X) %*% y)
    XbOLS    <- as.numeric(X %*% betaOLS)
    vOLS     <- var(y - XbOLS)
    Vinv_tmp <- solve(vOLS * h2 * K + vOLS * (1 - h2) * diag(n))
    betaGLS  <- solve(t(X) %*% Vinv_tmp %*% X) %*% (t(X) %*% Vinv_tmp %*% y)
    XbGLS    <- as.numeric(X %*% betaGLS)
    vGLS     <- var(y - XbGLS)
    sig2g    <- vGLS * h2
    sig2e    <- vGLS * (1 - h2)
  }

  V  <- sig2g * K + sig2e * diag(n)
  Vi <- solve(V)

  # Fixed-effect LOO estimates
  xtvxi    <- solve(t(X) %*% Vi %*% X) %*% (t(X) %*% Vi)
  betaFe   <- xtvxi %*% y
  Xproj    <- X %*% xtvxi
  yfe      <- Xproj %*% y
  yfeLoo   <- (yfe - diag(Xproj) * y) / (1 - diag(Xproj))

  yTilde   <- y - yfeLoo     # LOO residuals
  yFeResid <- y - yfe        # in-sample fixed-effect residuals

  # Intercept-only projection (used to re-centre in LOO folds)
  X1    <- matrix(1, n, 1)
  X1proj <- X1 %*% solve(t(X1) %*% Vi %*% X1) %*% (t(X1) %*% Vi)
  P1    <- diag(n) - X1proj

  # BLUPs
  bBlupStd <- (1 / m) * sig2g * t(Z) %*% (Vi %*% yFeResid)
  bBlup    <- (1 / m) * sig2g * t(Z) %*% (Vi %*% P1 %*% yTilde)
  H        <- sig2g * K %*% Vi %*% P1
  gBlup    <- H %*% yTilde
  gBlupLoo <- (gBlup - diag(H) * yTilde) / (1 - diag(H))

  ret <- list(
    betaFe        = betaFe,
    yfe           = yfe,
    yfeLoo        = yfeLoo,
    bBlup         = bBlup,
    bBlupStandard = bBlupStd,
    gBlup         = gBlup,
    gBlupLoo      = gBlupLoo
  )

  if (return_bBlupsLoo) {
    ret$bBlupsLoo <- .getBBlupsLoo(y = y, Z = Z, Vi = Vi,
                                   yfeLoo = yfeLoo,
                                   gBlupLoo = gBlupLoo,
                                   sig2g = sig2g)
  }

  if (return_intermediates) {
    ret$intermediates <- list(V = V, Vi = Vi, H = H,
                              sig2g = sig2g, sig2e = sig2e)
  }

  ret
}
