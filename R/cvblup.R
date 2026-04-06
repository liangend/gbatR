#' Compute cross-validated BLUPs for all genes on one chromosome
#'
#' For each gene on a given chromosome, extracts cis SNPs, builds a kinship
#' matrix, estimates heritability via REML, and returns leave-one-out BLUP
#' predictions (i.e. the cis-genetic component of expression for each gene).
#'
#' @param chr Chromosome identifier (integer or character, e.g. `1` or `"1"`).
#'   Used to substitute the `{chr}` placeholder in `bfile_pattern` and to
#'   filter rows in `gene_pos`.
#' @param gene_pos Data frame of gene annotations. Must contain a chromosome
#'   column named `"chr"` (values like `"chr1"`) and a position column named
#'   `"start"`. The gene ID column is specified by `gene_col`.
#' @param genotype_dir Path to the directory containing PLINK bfiles.
#' @param expr Normalized expression matrix, samples x genes. Column names
#'   must match the gene IDs in `gene_pos[[gene_col]]`.
#' @param gene_col Name of the gene ID column in `gene_pos` (default `"gene"`).
#' @param bfile_pattern Glob-style template for the PLINK bfile prefix within
#'   `genotype_dir`. Use `{chr}` as a placeholder for the chromosome value
#'   supplied in `chr`. Examples:
#'   \itemize{
#'     \item `"chr{chr}_QCed"` (default) matches `chr1_QCed`
#'     \item `"chr{chr}*"` matches any file starting with `chr1`
#'     \item `"*chr{chr}"` matches any file ending with `chr1`
#'     \item `"*chr{chr}*"` matches any file containing `chr1`
#'   }
#'   If the pattern matches multiple files the first is used.
#' @param svs Surrogate variable matrix (samples x K) from [compute_svs()],
#'   or `NA` if none.
#' @param covariates Additional covariate matrix (samples x C), or `NULL`.
#'   Combined with `svs` as fixed effects.
#' @param plink_path Path to PLINK 1.9 executable (default `"plink"`).
#' @param plink_samples Path to a two-column (FID IID) plain-text file
#'   specifying the samples to keep from the PLINK bfiles. If `NULL`, all
#'   samples in the bfile are used.
#' @param work_dir Directory for temporary PLINK output files
#'   (default `tempdir()`).
#' @param cis_window Window in base pairs around TSS used to define cis SNPs
#'   (default 100000).
#'
#' @return A list with:
#'   \item{blup}{N x G matrix of LOO BLUP predictions (G = genes processed).}
#'   \item{h2g}{Data frame with columns matching `gene_col` and `h2g`.}
#'
#' @export
run_cvblup <- function(chr, gene_pos, genotype_dir, expr,
                       gene_col       = "gene",
                       bfile_pattern  = "chr{chr}_QCed",
                       svs            = NA,
                       covariates     = NULL,
                       plink_path     = "plink",
                       plink_samples  = NULL,
                       work_dir       = tempdir(),
                       cis_window     = 1e5) {

  chr        <- as.character(chr)
  chrom_str  <- paste0("chr", chr)

  # ---- Resolve PLINK bfile via pattern ------------------------------------
  pattern_resolved <- gsub("{chr}", chr, bfile_pattern, fixed = TRUE)
  bim_matches <- Sys.glob(
    file.path(genotype_dir, paste0(pattern_resolved, ".bim"))
  )
  if (length(bim_matches) == 0) {
    stop("No PLINK bfile found matching pattern '", pattern_resolved,
         "' in directory: ", genotype_dir)
  }
  if (length(bim_matches) > 1) {
    warning("Multiple bfiles match pattern '", pattern_resolved,
            "'; using: ", bim_matches[1])
  }
  bfile    <- sub("\\.bim$", "", bim_matches[1])
  bim_file <- bim_matches[1]

  allmap     <- utils::read.table(bim_file, stringsAsFactors = FALSE)
  gene_names <- colnames(expr)

  # ---- All genes on this chromosome ---------------------------------------
  chr_flag <- which(gene_pos$chr == chrom_str)
  if (length(chr_flag) == 0) {
    message("No genes found for chromosome ", chr, " in gene_pos.")
    return(list(blup = NULL,
                h2g  = stats::setNames(
                  data.frame(character(0), numeric(0)),
                  c(gene_col, "h2g")
                )))
  }

  # ---- Build combined covariate matrix (SVs + covariates), standardised ---
  if (all(is.na(svs))) {
    sva_mat <- covariates
  } else {
    sv_mat  <- as.matrix(svs)
    sva_mat <- if (!is.null(covariates)) cbind(covariates, sv_mat) else sv_mat
  }

  if (!is.null(sva_mat)) {
    mean_covar <- colMeans(sva_mat)
    sva_mat    <- sweep(sva_mat, 2, mean_covar, "-")
    sd_covar   <- apply(sva_mat, 2, stats::sd)
    sd_covar[sd_covar == 0] <- 1
    sva_mat    <- sweep(sva_mat, 2, sd_covar, "/")
  }

  # ---- Compute sample alignment: PLINK order → expression order -----------
  fam_data   <- utils::read.table(paste0(bfile, ".fam"),
                                  stringsAsFactors = FALSE)
  plink_iids <- fam_data[, 2]
  if (!is.null(plink_samples)) {
    keep_iids  <- utils::read.table(plink_samples,
                                    stringsAsFactors = FALSE)[, 2]
    plink_iids <- plink_iids[plink_iids %in% keep_iids]
  }
  expr_ids  <- rownames(expr)
  align_idx <- match(expr_ids, plink_iids)
  if (anyNA(align_idx)) {
    stop("Some rownames(expr) not found in PLINK sample IDs. ",
         "Ensure rownames(expr) match the IID column of the .fam file.")
  }

  # ---- Per-gene processing ------------------------------------------------
  goos     <- NULL
  colgenes <- NULL
  h2_all   <- NULL

  prefix <- file.path(work_dir, paste0("chr", chr))

  for (flag_i in chr_flag) {
    gene       <- gene_pos[[gene_col]][flag_i]
    exp_gflag  <- which(gene_names == gene)
    if (length(exp_gflag) == 0) next

    tss1     <- max(0, gene_pos$start[flag_i] - cis_window)
    tss2     <- gene_pos$start[flag_i] + cis_window
    cis_flag <- which(allmap[, 4] > tss1 & allmap[, 4] < tss2)
    if (length(cis_flag) == 0) next

    # Write SNP range file for PLINK --extract --range
    snp_range_file <- paste0(prefix, "_snps.txt")
    utils::write.table(allmap[cis_flag, c(1, 4, 4, 4)],
                       snp_range_file,
                       quote = FALSE, row.names = FALSE, col.names = FALSE)

    # Run PLINK (quote all paths to handle spaces in directory names)
    plink_cmd <- paste(
      shQuote(plink_path),
      "--bfile", shQuote(bfile),
      "--extract", shQuote(snp_range_file), "--range",
      "--recode12",
      "--out", shQuote(prefix)
    )
    if (!is.null(plink_samples)) {
      plink_cmd <- paste(plink_cmd, "--keep", shQuote(plink_samples))
    }
    system(plink_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

    f_ped <- paste0(prefix, ".ped")
    f_map <- paste0(prefix, ".map")
    if (!file.exists(f_ped) || !file.exists(f_map)) next

    plink_out <- .read_plink_ped(f_ped, f_map)
    W <- .normalize_genotypes(plink_out$W[align_idx, , drop = FALSE])
    M <- ncol(W)
    if (M == 0) next

    K <- W %*% t(W) / M
    phe <- expr[, exp_gflag]

    pError <- tryCatch(
      aiML(A = list(K), y = phe, Var = c(0.5, 0.5),
           verbose = FALSE, CALC_SE = FALSE),
      error = function(e) e
    )

    if (!inherits(pError, "error")) {
      reml_est <- aiML(A = list(K), y = phe, Var = c(0.5, 0.5),
                       verbose = FALSE, CALC_SE = FALSE)
      h2g <- as.numeric(reml_est$h2)
    } else {
      h2g <- 0.05
    }

    blup_out <- .getBlups(y = phe, X = sva_mat, Z = W,
                          sig2g = NA, sig2e = NA, h2 = h2g)

    colgenes <- c(colgenes, gene)
    h2_all   <- c(h2_all, h2g)
    goos     <- cbind(goos, blup_out$gBlupLoo)
  }

  # ---- Clean up PLINK temp files ------------------------------------------
  for (ext in c(".ped", ".map", ".nosex", ".log", ".nof")) {
    f <- paste0(prefix, ext)
    if (file.exists(f)) file.remove(f)
  }
  snp_range_file <- paste0(prefix, "_snps.txt")
  if (file.exists(snp_range_file)) file.remove(snp_range_file)

  if (is.null(goos)) {
    return(list(blup = NULL,
                h2g  = stats::setNames(
                  data.frame(character(0), numeric(0)),
                  c(gene_col, "h2g")
                )))
  }

  colnames(goos) <- colgenes
  h2g_df <- data.frame(colgenes, h2_all, stringsAsFactors = FALSE)
  colnames(h2g_df) <- c(gene_col, "h2g")

  list(blup = goos, h2g = h2g_df)
}
