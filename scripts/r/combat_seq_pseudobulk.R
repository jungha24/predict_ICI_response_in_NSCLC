args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: Rscript combat_seq_pseudobulk.R <counts_csv> <meta_csv> <output_csv> <batch_col> <covariates_csv>")
}

counts_csv <- args[[1]]
meta_csv <- args[[2]]
output_csv <- args[[3]]
batch_col <- args[[4]]
covariates_csv <- args[[5]]

suppressPackageStartupMessages({
  library(sva)
})

counts <- read.csv(counts_csv, row.names = 1, check.names = FALSE)
meta <- read.csv(meta_csv, check.names = FALSE, stringsAsFactors = FALSE)

if (!(batch_col %in% colnames(meta))) {
  stop(sprintf("batch column '%s' was not found in metadata.", batch_col))
}

counts_mat <- as.matrix(round(counts))
batch <- meta[[batch_col]]

covariates <- character(0)
if (nzchar(covariates_csv)) {
  covariates <- trimws(unlist(strsplit(covariates_csv, ",")))
  covariates <- covariates[nzchar(covariates)]
}

covar_mod <- NULL
if (length(covariates) > 0) {
  missing_covars <- setdiff(covariates, colnames(meta))
  if (length(missing_covars) > 0) {
    stop(sprintf("covariates missing from metadata: %s", paste(missing_covars, collapse = ", ")))
  }

  covar_mod <- meta[, covariates, drop = FALSE]
  for (col in colnames(covar_mod)) {
    x <- covar_mod[[col]]
    numeric_x <- suppressWarnings(as.numeric(x))
    if (all(is.na(x) == is.na(numeric_x))) {
      covar_mod[[col]] <- numeric_x
      med <- median(covar_mod[[col]], na.rm = TRUE)
      if (is.na(med)) med <- 0
      covar_mod[[col]][is.na(covar_mod[[col]])] <- med
    } else {
      x[is.na(x) | x == ""] <- "NA"
      covar_mod[[col]] <- factor(x)
    }
  }
}

adjusted <- ComBat_seq(
  counts = counts_mat,
  batch = batch,
  covar_mod = covar_mod
)

write.csv(adjusted, output_csv, quote = FALSE)
