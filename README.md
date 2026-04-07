# gbatR

**Gene-Based Association Test for Trans-eQTL Discovery**

gbatR identifies trans-associations between genes by modelling the cis-genetic component of each gene's expression as a predictor and testing it against all trans-genes. Please check Liu et al. (2020) for details: <https://doi.org/10.1186/s13059-020-02120-1>

------------------------------------------------------------------------

## Table of Contents

1.  [Dependencies](#dependencies)
2.  [Installation](#installation)
3.  [Input Data Requirements](#input-data-requirements)
4.  [Run the Full Pipeline in One Call](#run-the-full-pipeline-in-one-call)
5.  [Run the Pipeline Step by Step](#run-the-pipeline-step-by-step)
6.  [Parallel Chromosome Processing](#parallel-chromosome-processing)
7.  [Output Files](#output-files)

------------------------------------------------------------------------

## Dependencies

### R packages

| Package | Source | Role |
|------------------------|------------------------|------------------------|
| `limma` | Bioconductor | Quantile normalization |
| `sva` | Bioconductor | Surrogate variable analysis |
| `qvalue` | Bioconductor | FDR estimation |
| `RNOmni` | CRAN | Rank-based inverse-normal transform |
| `msm` | CRAN | Delta method for REML standard errors |
| `SmartSVA` | CRAN | Fast SVA for trans p-value computation |
| `phenix` | <https://mathgen.stats.ox.ac.uk/genetics_software/phenix/phenix.html> | Quantile normalization of residuals |

Install all dependencies:

``` r
# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("limma", "sva", "qvalue"))

# CRAN packages
install.packages(c("RNOmni", "msm", "SmartSVA"))

# phenix package. Please download the package from the link provided above to a local path
install.packages(path_to_phenix, repos = NULL, type="source")
```

### External software

**PLINK** is required for extracting cis-SNP genotypes. Download from <https://www.cog-genomics.org/plink/1.9/> and note the path to the executable.

------------------------------------------------------------------------

## Installation

``` r
remotes::install_github("liangend/gbatR")
library(gbat)
```

------------------------------------------------------------------------

## Input Data Requirements

| Object | Format | Notes |
|------------------------|------------------------|------------------------|
| `expr` | Gene expression matrix, samples × genes | Row names are sample IDs and column names are genes. |
| `gene_info` | Gene information data frame | Must contain a gene name column (matching with gene names in expr), a `chr` column in `"chr1"` format, and a `start` column (TSS position in the same coordinate system as the PLINK bfiles). |
| `covariates` | Covariate matrix, samples × covariates | Same row order as `expr`. Typical columns: genotype PCs, age, sex, batch. Pass `NULL` if none. |
| PLINK bfiles | `.bed` / `.bim` / `.fam` | One set per chromosome. File names must follow a consistent pattern with `{chr}` as the chromosome placeholder, e.g. ``` chr1_eur.bed,``chr2_eur.bed ```. |
| `plink_samples` | Plain-text file, two columns (FID IID, no header) | Samples overlapping between gene expression and genotype. |

> **Sample alignment:** `rownames(expr)` must match the `rownames(covariates)`exactly.

> **Coordinate system:** The `start` column in `gene_info` and the SNP positions in the PLINK `.bim` file must use the same genome build (e.g. both hg38 or both hg19).

------------------------------------------------------------------------

## Run the Full Pipeline in One Call

`run_gbat()` executes steps 1–5 (normalization through p-value aggregation) and writes all intermediate results to `output_dir`. FDR computation and significant-gene extraction are then run separately (see below). In the `example` directory, we provided a toy example to demonstrate the pipeline.

`run_gbat()` accepts either normalized gene expression matrix or raw gene count matrix. When `expr` is a gene count matrix, `run_gbat()` can normalize the the matrix by setting `normalize = TRUE`. Under this case, `gene_info` must contain a column of `length` recording the gene length (the total number of nucleotides covered by its exons), which can be obtained through gene count software (e.g., featureCounts, htseq-count). The `length` information is required to calculate RPKM in the normalization procedure.

``` r
library(gbat)

# Load inputs
expr      <- readRDS("example/gene_expr.rds")                          # samples × genes
gene_info <- read.table("example/gene_info.txt", header = T)           # gene annotation table
cov       <- read.table("example/covariates.txt",
                        header = TRUE, row.names = 1, sep = "\t")
output_dir <- "~/gbat_output"  # customize it to whatever you want

# Extract the overlapping IDs between expr and genotype (IMPORTANT!)
plink_samples <- read.table("example/chr1_QCed.fam")
keep_ids <- intersect(plink_samples$V2, rownames(expr))
keep_file <- file.path(output_dir, "plink_keep.txt")
FID <- plink_samples$V1[match(keep_ids, plink_samples$V2)]
write.table(data.frame(FID = FID, IID = keep_ids),
            keep_file, quote = FALSE, row.names = FALSE, col.names = FALSE)

# Keep overlapping IDs
expr <- expr[keep_ids, ]
cov <- cov[keep_ids, ]

# Run pipeline (steps 1–5)
results <- run_gbat(
  expr          = expr,
  gene_info     = gene_info,
  covariates    = as.matrix(cov),
  genotype_dir  = "example/",
  output_dir    = output_dir,
  gene_col      = "gene_id",          # name of gene name column in gene_info that matches with expr
  bfile_pattern = "chr{chr}_QCed",    # PLINK file name format. {chr} is the chromosome number
  plink_path    = "/path/to/plink",
  plink_samples = keep_file,          # NULL if PLINK file contains the same IDs as expr. The ID order doesn't have to match with expr
  num_sv        = 3,                  # number of surrogated variables. default = 20
  cis_window    = 1e5,                # cis window around TSS to predict the gene expression
  chromosomes   = 1:5,                # better to run one chr at a time in real data analysis
  normalize     = FALSE               # set TRUE if expr is a raw count matrix
)
```

**Step 6 – Genome-wide FDR:**

Genome-wide FDR cutoff is calculated based on the p values of inter-chromosmome gene pairs.

``` r
allp_files <- list.files("gbat_output",
                         pattern = "_allp_h2g_filtered_interchrom\\.txt$",
                         full.names = TRUE)
pval_list  <- lapply(allp_files, function(f) as.numeric(read.table(f)[, 1]))
## two FDR cutoffs
cutoffs    <- compute_qvals(pval_list, fdr_levels = c(0.05, 0.1))
print(cutoffs)
```

**Step 7 – Significant trans-eQTL gene pairs:**

``` r
chrs <- 1:5

cor_list <- setNames(lapply(chrs, function(chr) {
  f <- file.path("gbat_output", paste0("cor_chr", chr, ".txt"))
  if (file.exists(f)) read.table(f, header = TRUE, row.names = 1) else NULL
}), as.character(chrs))

sig <- get_sig_genes(
  cor_results_list = cor_list,
  gene_pos         = gene_info,
  cutoff           = cutoffs["fdr0.05"],
  gene_col         = "gene_id"
)

write.table(sig, "gbat_output/sig_trans_eqtl.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)
head(sig)
```

------------------------------------------------------------------------

## Run the Pipeline Step by Step

The same pipeline can be assembled from individual functions, which is useful for custom workflows.

### Step 0 (optional) — Normalize raw counts

Skip this step if `expr` is already normalized.

``` r
norm_out  <- normalize_counts(count_matrix, gene_info,
                              min_cpm = 0.5, min_sample_frac = 0.5)
expr      <- norm_out$expr       # samples × genes, quantile-normalized
gene_info <- norm_out$gene_info  # filtered to expressed autosomal genes
```

`normalize_counts()` removes sex/mitochondrial genes, filters lowly expressed genes (CPM threshold), computes RPKM, and applies two rounds of quantile normalization.

### Step 1 — Surrogate variable analysis

``` r
svs <- compute_svs(expr, num_sv = 20)
# Returns a samples × num_sv matrix.
# Save to disk so parallel chromosome jobs can reuse it:
write.table(svs, "gbat_output/sva.txt", quote = FALSE, row.names = FALSE)
```

### Step 2 — Cross-validated BLUPs (per chromosome)

For each gene on the target chromosome, `run_cvblup()` extracts cis-SNPs, builds a kinship matrix, estimates heritability via REML, and returns leave-one-out BLUP predictions (the cis-genetic component of expression).

``` r
blup_out <- run_cvblup(
  chr           = 1,
  gene_pos      = gene_info,
  genotype_dir  = "/path/to/plink_bfiles",
  expr          = expr,
  gene_col      = "gene_id",
  bfile_pattern = "chr{chr}_QCed",
  svs           = svs,
  covariates    = as.matrix(cov),
  plink_path    = "/path/to/plink",
  plink_samples = "plink_keep.txt",
  work_dir      = "gbat_output",
  cis_window    = 1e5
)

# blup_out$blup  : samples × genes BLUP matrix
# blup_out$h2g   : data frame of per-gene heritability estimates
```

### Step 3 — Pearson R quality control

Measures how well each gene's BLUP predicts its own expression. Genes with Pearson R ≤ 0.1 are skipped in the trans p-value step.

``` r
pr <- compute_pearson_r(blup_out$blup, expr)
# Returns a data frame with columns: gene, pearsonR
```

### Step 4 — Trans p-values

For each gene with Pearson R \> 0.1, re-estimates surrogate variables via SmartSVA (using the BLUP as a covariate), residualizes expression, and regresses all genes on the standardized BLUP.

``` r
pvals <- compute_trans_pvals(expr, blup_out$blup, pr, num_sv = 20)
# Returns a genes × predictor-genes p-value matrix.
# Row names are target gene IDs; column names are predictor gene IDs.

write.table(pvals, "gbat_output/cor_chr1.txt",
            row.names = TRUE, col.names = TRUE, quote = FALSE)
```

### Step 5 — Aggregate inter-chromosomal p-values

Removes same-chromosome associations and genes with invalid heritability estimates (h2g outside [0, 1]), returning a flat vector of trans p-values for the predictor genes on this chromosome.

``` r
allp_chr1 <- aggregate_pvals(
  cor_results_list = list(pvals),
  h2g_results_list = list(blup_out$h2g),
  gene_pos         = gene_info,
  chr              = 1,
  gene_col         = "gene_id"
)

write.table(allp_chr1,
            "gbat_output/chr1_allp_h2g_filtered_interchrom.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
```

### Step 6 — Genome-wide FDR

Repeat steps 2–5 for all chromosomes, then pool all aggregated p-values:

``` r
allp_files <- list.files("gbat_output",
                         pattern = "_allp_h2g_filtered_interchrom\\.txt$",
                         full.names = TRUE)
pval_list  <- lapply(allp_files, function(f) as.numeric(read.table(f)[, 1]))
cutoffs    <- compute_qvals(pval_list, fdr_levels = c(0.05, 0.1))
print(cutoffs)
```

### Step 7 — Significant trans-eQTL gene pairs

``` r
# Read per-chromosome p-value matrices from disk (written in step 4)
cor_list <- setNames(lapply(1:22, function(chr) {
  f <- file.path("gbat_output", paste0("cor_chr", chr, ".txt"))
  if (file.exists(f)) read.table(f, header = TRUE, row.names = 1) else NULL
}), as.character(1:22))

sig <- get_sig_genes(
  cor_results_list     = cor_list,
  gene_pos             = gene_info,
  cutoff               = cutoffs["fdr0.05"],
  gene_col             = "gene_id",
  cis_exclusion_window = 1e7          # exclude pairs within 10 Mb on same chr
)
```

------------------------------------------------------------------------

## Output Files

All files are written to `output_dir`.

| File | Description |
|------------------------------------|------------------------------------|
| `sva.txt` | Surrogate variable matrix (samples × num_sv) |
| `chr{N}_h2g.txt` | Per-gene heritability estimates for chromosome N |
| `pearsonR_chr{N}.txt` | Pearson R between BLUP and observed expression for chromosome N |
| `cor_chr{N}.txt` | Trans p-value matrix for chromosome N (genes × predictor genes); row names are target gene IDs |
| `chr{N}_allp_h2g_filtered_interchrom.txt` | Flat vector of inter-chromosomal trans p-values for chromosome N (input to FDR) |
| `qval_cutoffs.txt` | P-value cutoffs at each requested FDR level |
| `sig_trans_eqtl.txt` | Final table of significant trans-eQTL gene pairs |

### Significant trans-eQTL output columns

One row per significant regulator–target pair, ordered by regulator chromosome and start position.

| Column | Description |
|------------------------------------|------------------------------------|
| `regulator_gene` | Regulator gene ID (the gene whose cis-genetic component predicts the target) |
| `regulator_chr` | Regulator gene chromosome |
| `regulator_start` / `regulator_end` | Regulator gene genomic coordinates |
| `target_gene` | Target gene ID (the trans-regulated gene) |
| `target_chr` | Target gene chromosome |
| `target_start` / `target_end` | Target gene genomic coordinates |
| `pval` | GBAT trans p-value for this regulator–target pair |
