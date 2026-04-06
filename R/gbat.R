#' Run the full GBAT trans-eQTL discovery pipeline
#'
#' Starting from a (pre-)normalized expression matrix, this function:
#' 1. Optionally normalizes expression from raw counts (RPKM + double
#'    quantile normalization). Skip by passing an already-normalized matrix
#'    directly as `expr`.
#' 2. Estimates surrogate variables (SVA).
#' 3. For each chromosome, computes cross-validated BLUPs and Pearson R
#'    diagnostics.
#' 4. Computes trans p-values via SmartSVA residualization.
#' 5. Aggregates inter-chromosomal p-values.
#' 6. Computes genome-wide FDR thresholds with `qvalue`.
#' 7. Returns significant trans-eQTL gene pairs.
#'
#' @param expr Normalized expression matrix (samples x genes), **or** a raw
#'   count matrix (genes x samples) when `normalize = TRUE`.
#' @param gene_info Data frame of gene annotations. Must contain columns for
#'   chromosome (`"chr"`), TSS position (`"start"`), and gene ID (specified
#'   by `gene_col`). Chromosome values must be in `"chr1"` format.
#' @param covariates Numeric matrix of additional covariates (samples x C),
#'   e.g. genotype PCs, in the same sample order as `expr`. Pass `NULL` if
#'   none.
#' @param genotype_dir Path to the directory containing PLINK bfiles.
#' @param output_dir Directory where intermediate and final results are
#'   written. Created if it does not exist.
#' @param gene_col Name of the gene ID column in `gene_info` (default
#'   `"gene"`). Column values must match column names of `expr`.
#' @param bfile_pattern Glob-style template for the PLINK bfile prefix
#'   within `genotype_dir`. Use `{chr}` as a placeholder for the chromosome.
#'   Default `"chr{chr}_QCed"`. See [run_cvblup()] for examples.
#' @param plink_path Path to PLINK 1.9 executable (default `"plink"`).
#' @param plink_samples Path to a PLINK `--keep` file (FID IID, no header).
#'   If `NULL`, all samples in the bfile are used.
#' @param num_sv Number of surrogate variables (default 20).
#' @param cis_window Cis window around TSS in bp (default 100000).
#' @param chromosomes Integer vector of chromosomes to analyse (default 1:22).
#' @param normalize If `TRUE`, treat `expr` as a raw count matrix and run
#'   [normalize_counts()] first. Requires `gene_info` to have a `length`
#'   column. Default `FALSE`.
#' @param min_cpm Minimum CPM for expression filter (only used when
#'   `normalize = TRUE`, default 0.5).
#' @param min_sample_frac Minimum sample fraction for expression filter (only
#'   used when `normalize = TRUE`, default 0.5).
#' @param pseudogenes Optional character vector of gene IDs to exclude.
#' @param verbose Print progress messages (default TRUE).
#'
#' @return A list with:
#'   \item{expr}{Expression matrix used (samples x genes).}
#'   \item{gene_info}{Gene annotation (filtered if `normalize = TRUE`).}
#'   \item{svs}{Surrogate variable matrix.}
#'   \item{blup_results}{List of BLUP outputs, one entry per chromosome.}
#'
#' @export
run_gbat <- function(expr,
                     gene_info,
                     covariates      = NULL,
                     genotype_dir,
                     output_dir,
                     gene_col        = "gene",
                     bfile_pattern   = "chr{chr}_QCed",
                     plink_path      = "plink",
                     plink_samples   = NULL,
                     num_sv          = 20,
                     cis_window      = 1e5,
                     chromosomes     = 1:22,
                     normalize       = FALSE,
                     min_cpm         = 0.5,
                     min_sample_frac = 0.5,
                     pseudogenes     = NULL,
                     verbose         = TRUE) {

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # ---- Step 1: Normalize expression (optional) -----------------------------
  if (normalize) {
    if (verbose) message("[GBAT] Step 1: Normalizing expression...")
    norm_out  <- normalize_counts(expr, gene_info,
                                  min_cpm = min_cpm,
                                  min_sample_frac = min_sample_frac)
    expr      <- norm_out$expr
    gene_info <- norm_out$gene_info
    saveRDS(list(expr = expr, gene_info = gene_info),
            file.path(output_dir, "expr_normalized.rds"))
  } else {
    if (verbose) message("[GBAT] Step 1: Using pre-normalized expression matrix.")
  }

  # ---- Step 2: Surrogate variable analysis ---------------------------------
  sva_file <- file.path(output_dir, "sva.txt")
  if (file.exists(sva_file)) {
    if (verbose) message("[GBAT] Step 2: Loading existing surrogate variables from ", sva_file)
    svs <- as.matrix(utils::read.table(sva_file, header = TRUE))
  } else {
    if (verbose) message("[GBAT] Step 2: Computing surrogate variables...")
    svs <- compute_svs(expr, num_sv = num_sv)
    utils::write.table(svs, sva_file, quote = FALSE, row.names = FALSE)
  }

  # ---- Steps 3–5: Per-chromosome processing --------------------------------
  blup_all     <- list()
  pearsonr_all <- list()
  pval_all     <- list()
  h2g_all      <- list()

  for (chr in chromosomes) {
    chr_str <- as.character(chr)
    if (verbose) message("[GBAT] Chr ", chr, ": cvBLUP + trans p-values...")

    # Check bfile exists before running
    pattern_resolved <- gsub("{chr}", chr_str, bfile_pattern, fixed = TRUE)
    bim_check <- Sys.glob(file.path(genotype_dir,
                                    paste0(pattern_resolved, ".bim")))
    if (length(bim_check) == 0) {
      if (verbose) message("  Skipping chr", chr, ": no bfile matching '",
                           pattern_resolved, "' found.")
      next
    }

    blup_out <- run_cvblup(
      chr           = chr,
      gene_pos      = gene_info,
      genotype_dir  = genotype_dir,
      expr          = expr,
      gene_col      = gene_col,
      bfile_pattern = bfile_pattern,
      svs           = svs,
      covariates    = covariates,
      plink_path    = plink_path,
      plink_samples = plink_samples,
      work_dir      = output_dir,
      cis_window    = cis_window
    )

    if (is.null(blup_out$blup) || ncol(blup_out$blup) == 0) {
      if (verbose) message("  Chr ", chr, ": no genes with valid BLUPs — ",
                           "check that cis-SNPs exist within +/-", cis_window,
                           " bp of each gene's TSS and that sample IDs match.")
      blup_all[[chr_str]]     <- NULL
      pearsonr_all[[chr_str]] <- NULL
      pval_all[[chr_str]]     <- NULL
      h2g_all[[chr_str]]      <- blup_out$h2g
      next
    }
    if (verbose) message("  Chr ", chr, ": BLUPs computed for ",
                         ncol(blup_out$blup), " genes.")

    utils::write.table(blup_out$h2g,
                       file.path(output_dir, paste0("chr", chr, "_h2g.txt")),
                       col.names = TRUE, row.names = FALSE,
                       sep = "\t", quote = FALSE)

    pr <- compute_pearson_r(blup_out$blup, expr)
    n_pass_pr <- sum(pr$pearsonR > 0.1, na.rm = TRUE)
    if (verbose) message("  Chr ", chr, ": ", n_pass_pr, " / ",
                         nrow(pr), " genes pass pearsonR > 0.1 threshold.")
    if (n_pass_pr == 0) {
      if (verbose) message("  Chr ", chr, ": skipping trans p-values — ",
                           "no genes with sufficient cis-genetic signal ",
                           "(pearsonR > 0.1). Consider whether h2g is ",
                           "estimable with the current sample size and ",
                           "cis-window.")
      utils::write.table(pr,
                         file.path(output_dir, paste0("pearsonR_chr", chr, ".txt")),
                         row.names = FALSE, quote = FALSE, sep = "\t")
      blup_all[[chr_str]]     <- blup_out$blup
      pearsonr_all[[chr_str]] <- pr
      pval_all[[chr_str]]     <- NULL
      h2g_all[[chr_str]]      <- blup_out$h2g
      next
    }
    utils::write.table(pr,
                       file.path(output_dir, paste0("pearsonR_chr", chr, ".txt")),
                       row.names = FALSE, quote = FALSE, sep = "\t")

    pvals <- compute_trans_pvals(expr, blup_out$blup, pr, num_sv = num_sv)
    if (ncol(pvals) == 0) {
      if (verbose) message("  Chr ", chr, ": trans p-value computation ",
                           "produced no results — SmartSVA may have failed ",
                           "for all genes. Check that num_sv (", num_sv,
                           ") is appropriate for the number of genes (",
                           ncol(expr), ").")
      blup_all[[chr_str]]     <- blup_out$blup
      pearsonr_all[[chr_str]] <- pr
      pval_all[[chr_str]]     <- NULL
      h2g_all[[chr_str]]      <- blup_out$h2g
      next
    }
    utils::write.table(pvals,
                       file.path(output_dir, paste0("cor_chr", chr, ".txt")),
                       row.names = TRUE, col.names = colnames(pvals))

    blup_all[[chr_str]]     <- blup_out$blup
    pearsonr_all[[chr_str]] <- pr
    pval_all[[chr_str]]     <- pvals
    h2g_all[[chr_str]]      <- blup_out$h2g

    chr_allp <- aggregate_pvals(
      list(pvals), list(blup_out$h2g), gene_info, chr,
      gene_col = gene_col
    )
    if (length(chr_allp) == 0) {
      if (verbose) message("  Chr ", chr, ": no inter-chromosomal p-values ",
                           "after filtering — check that gene_info contains ",
                           "genes on chromosomes other than chr", chr, " and ",
                           "that gene IDs in gene_info match rownames(expr).")
    } else {
      if (verbose) message("  Chr ", chr, ": ", length(chr_allp),
                           " inter-chromosomal p-values aggregated.")
    }
    utils::write.table(
      chr_allp,
      file.path(output_dir,
                paste0("chr", chr, "_allp_h2g_filtered_interchrom.txt")),
      quote = FALSE, row.names = FALSE, col.names = FALSE
    )
  }

  list(
    expr         = expr,
    gene_info    = gene_info,
    svs          = svs,
    blup_results = blup_all
  )
}
