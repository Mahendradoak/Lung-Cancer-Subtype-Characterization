# Multi-Omics Factor Analysis. The original repository only contained a script
# that READ a pre-trained model ("trained_mofa_model.rds"), which was never
# committed, so the factors could not be regenerated. This trains it.

top_variable <- function(m, n) {
  v <- matrixStats::rowVars(m, na.rm = TRUE)
  m[order(v, decreasing = TRUE)[seq_len(min(n, nrow(m)))], , drop = FALSE]
}

train_mofa <- function(expr, cnv, out_dir, n_expr = 5000, n_cnv = 5000, n_factors = 15, seed = 42) {
  if (!requireNamespace("MOFA2", quietly = TRUE)) {
    stop("MOFA2 is required to train factors: BiocManager::install('MOFA2'). ",
         "Or pass --factors with an existing factor matrix.")
  }
  views <- list(expression = top_variable(expr, n_expr))
  if (!is.null(cnv)) views$cnv <- top_variable(cnv, n_cnv)
  views <- lapply(views, function(m) m - rowMeans(m, na.rm = TRUE))
  obj <- MOFA2::create_mofa(views)
  dopts <- MOFA2::get_default_data_options(obj); dopts$scale_views <- TRUE
  mopts <- MOFA2::get_default_model_options(obj); mopts$num_factors <- n_factors
  topts <- MOFA2::get_default_training_options(obj); topts$seed <- seed; topts$convergence_mode <- "slow"
  obj <- MOFA2::prepare_mofa(obj, data_options = dopts, model_options = mopts, training_options = topts)
  model_file <- file.path(out_dir, "mofa_model.hdf5")
  obj <- MOFA2::run_mofa(obj, outfile = model_file, use_basilisk = TRUE)
  factors <- MOFA2::get_factors(obj, factors = "all")[[1]]
  r2 <- MOFA2::get_variance_explained(obj)$r2_per_factor[[1]]
  utils::write.csv(r2, file.path(out_dir, "mofa_variance_explained.csv"))
  factors
}

# How much of each factor is explained by LUAD vs LUSC histology (R^2 of a
# one-way linear model). Factors that mostly encode histology make any
# clustering on them rediscover histology.
factor_histology_r2 <- function(factors, histology) {
  ok <- !is.na(histology)
  data.frame(factor = colnames(factors),
             r2 = apply(factors[ok, , drop = FALSE], 2, function(f) summary(stats::lm(f ~ histology[ok]))$r.squared))
}

select_factors <- function(factors, histology, max_r2 = NULL) {
  assoc <- factor_histology_r2(factors, histology)
  keep <- if (is.null(max_r2)) assoc$factor else assoc$factor[assoc$r2 <= max_r2]
  if (length(keep) < 2) stop("fewer than 2 factors left after removing histology-associated factors")
  list(factors = factors[, keep, drop = FALSE], association = assoc, kept = keep)
}
