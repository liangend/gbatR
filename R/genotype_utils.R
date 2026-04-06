# Genotype matrix utilities: reading PLINK ped/map files and normalizing.


# Read a PLINK --recode12 ped/map pair into a genotype matrix.
# Returns a list: W (N x M matrix), fam (sample info), map (SNP info).
.read_plink_ped <- function(f_ped, f_map) {
  map <- utils::read.table(f_map, stringsAsFactors = FALSE, fill = TRUE)
  ped <- utils::read.table(f_ped, stringsAsFactors = FALSE)

  N <- nrow(ped)
  M <- nrow(map)
  W_all <- matrix(0, nrow = N, ncol = M)
  fam <- ped[, 1:6]

  for (s in seq_len(M)) {
    pos <- c(6 + s * 2 - 1, 6 + s * 2)
    al <- unique(as.integer(c(ped[, pos[1]], ped[, pos[2]])))
    al <- al[al != 0]
    for (r in pos) {
      W_all[ped[, r] != al[1], s] <- W_all[ped[, r] != al[1], s] + 1
    }
    W_all[ped[, pos[2]] == 0, s] <- NA
  }

  list(W = W_all, fam = fam, map = map)
}


# Standardize a genotype matrix: center by 2*MAF, scale by SD, zero-out
# monomorphics and missing values.
# Returns the normalized matrix W (N x M).
.normalize_genotypes <- function(W) {
  N <- nrow(W)
  M <- ncol(W)

  mafs <- colMeans(W, na.rm = TRUE) / 2

  # Center
  W_new <- sweep(W, 2, 2 * mafs, "-")

  # Set missing to zero after centering
  W_new[is.na(W_new)] <- 0

  # Scale to unit variance
  vars <- apply(W_new, 2, var)
  W_new <- sweep(W_new, 2, sqrt(vars), "/")

  # Zero out monomorphic columns (var == 0 -> Inf/NaN after division)
  W_new[is.na(W_new)] <- 0

  W_new
}
