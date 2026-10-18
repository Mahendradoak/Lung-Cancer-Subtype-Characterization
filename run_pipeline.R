#!/usr/bin/env Rscript
# Multi-omics subtyping of TCGA lung cancer (LUAD + LUSC).
#
# Usage:
#   Rscript run_pipeline.R --expr HiSeqV2 --cnv Gistic2_all_data_by_genes \
#       --clinical TCGA.LUNG.sampleMap_LUNG_clinicalMatrix --out results
# Options:
#   --factors FILE        reuse a factor matrix (rows = samples) instead of training MOFA
#   --max-histology-r2 X  drop factors whose R^2 with LUAD/LUSC exceeds X (e.g. 0.3)
#   --k K                 force the number of clusters (default: chosen by consensus PAC)
#   --adjust-histology    adjust differential expression for histology
#   --seed N

if (!exists("harmonise", mode = "function")) {   # not already loaded (e.g. by the tests)
  script_dir <- local({
    f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    if (length(f)) dirname(normalizePath(f)) else getwd()
  })
  for (f in list.files(file.path(script_dir, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)
}

parse_args <- function(args) {
  out <- list(seed = 42, k = NA, max_histology_r2 = NA, adjust_histology = FALSE, cnv = NA, factors = NA)
  i <- 1
  while (i <= length(args)) {
    key <- gsub("-", "_", sub("^--", "", args[i]))
    if (key == "adjust_histology") { out$adjust_histology <- TRUE; i <- i + 1; next }
    out[[key]] <- args[i + 1]; i <- i + 2
  }
  for (req in c("expr", "clinical", "out")) if (is.null(out[[req]])) stop("missing --", req)
  out$seed <- as.integer(out$seed); out$k <- as.integer(out$k)
  out$max_histology_r2 <- as.numeric(out$max_histology_r2)
  out
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  a <- parse_args(args)
  dir.create(a$out, showWarnings = FALSE, recursive = TRUE)
  set.seed(a$seed)

  message("1/7 loading and harmonising")
  cnv <- if (!is.na(a$cnv)) read_matrix(a$cnv) else NULL
  if (is.null(cnv)) message("   no --cnv given: expression-only analysis (results will say so)")
  h <- harmonise(read_matrix(a$expr), cnv, read_clinical(a$clinical))
  message(sprintf("   %d primary tumours (%s)", ncol(h$expr), paste(names(table(h$clinical$histology)),
                                                                     table(h$clinical$histology), collapse = ", ")))

  message("2/7 latent factors")
  factors <- if (!is.na(a$factors)) {
    f <- read_matrix(a$factors)
    f[intersect(rownames(f), colnames(h$expr)), , drop = FALSE]
  } else {
    f <- train_mofa(h$expr, h$cnv, a$out, seed = a$seed)
    utils::write.csv(f, file.path(a$out, "mofa_factors.csv"))
    f
  }
  keep_samples <- intersect(colnames(h$expr), rownames(factors))
  h$expr <- h$expr[, keep_samples]; if (!is.null(h$cnv)) h$cnv <- h$cnv[, keep_samples]
  h$clinical <- h$clinical[match(keep_samples, h$clinical$sampleID), ]
  factors <- factors[keep_samples, , drop = FALSE]
  sel <- select_factors(factors, h$clinical$histology, if (is.na(a$max_histology_r2)) NULL else a$max_histology_r2)
  utils::write.csv(sel$association, file.path(a$out, "factor_histology_r2.csv"), row.names = FALSE)

  message("3/7 choosing k")
  ks <- choose_k(sel$factors, seed = a$seed, out_dir = a$out)
  k <- if (!is.na(a$k)) a$k else ks$k
  clusters <- paste0("S", cluster_kmeans(sel$factors, k, a$seed))
  stab <- bootstrap_stability(sel$factors, k, seed = a$seed)
  hist_ag <- histology_agreement(clusters, h$clinical$histology)
  assign_df <- data.frame(sampleID = keep_samples, subtype = clusters, histology = h$clinical$histology)
  utils::write.csv(assign_df, file.path(a$out, "subtype_assignments.csv"), row.names = FALSE)
  utils::write.csv(as.data.frame.matrix(hist_ag$table), file.path(a$out, "subtype_vs_histology.csv"))

  message("4/7 differential expression, CNV, enrichment, targets")
  hist_cov <- if (a$adjust_histology) h$clinical$histology else NULL
  de_summary <- list()
  for (s in sort(unique(clusters))) {
    de <- de_one_vs_rest(h$expr, clusters, s, hist_cov)
    utils::write.csv(de, file.path(a$out, sprintf("%s_DE_limma.csv", s)), row.names = FALSE)
    cnv_res <- NULL
    if (!is.null(h$cnv)) {
      cnv_res <- diff_cnv(h$cnv, clusters, s)
      utils::write.csv(cnv_res, file.path(a$out, sprintf("%s_CNV_wilcoxon.csv", s)), row.names = FALSE)
    }
    en <- enrich(de)
    for (nm in names(en)) utils::write.csv(en[[nm]], file.path(a$out, sprintf("%s_%s.csv", s, nm)), row.names = FALSE)
    utils::write.csv(annotate_targets(de, cnv_res), file.path(a$out, sprintf("%s_target_annotation.csv", s)),
                     row.names = FALSE)
    de_summary[[s]] <- c(n_tested = nrow(de), n_sig = nrow(significant(de)))
    grDevices::pdf(file.path(a$out, sprintf("%s_volcano.pdf", s)), 7, 6)
    with(de, plot(log2FC, -log10(adj_p_value), pch = 20, cex = 0.4,
                  col = ifelse(adj_p_value < 0.05 & abs(log2FC) >= 1, "firebrick", "grey60"),
                  main = sprintf("%s vs rest (limma%s)", s, if (a$adjust_histology) ", histology-adjusted" else "")))
    grDevices::dev.off()
  }
  if (!is.null(h$cnv)) {
    utils::write.csv(cnv_expression_correlation(h$expr, h$cnv), file.path(a$out, "cnv_expression_spearman.csv"),
                     row.names = FALSE)
  }

  message("5/7 survival")
  sv <- survival_analysis(h$clinical, clusters)
  utils::write.csv(sv$cox, file.path(a$out, "cox_model.csv"), row.names = FALSE)
  grDevices::pdf(file.path(a$out, "kaplan_meier.pdf"), 7, 6)
  fit <- survival::survfit(survival::Surv(os_months, os_event) ~ clusters, data = h$clinical)
  plot(fit, col = seq_len(k), xlab = "Months", ylab = "Overall survival", main = sprintf("log-rank p = %.3g", sv$logrank_p))
  graphics::legend("topright", legend = sort(unique(clusters)), col = seq_len(k), lty = 1)
  grDevices::dev.off()

  message("6/7 summary")
  summary <- list(
    n_samples = length(keep_samples), cnv_used = !is.null(h$cnv), factors_used = sel$kept,
    k = k, k_chosen_by = if (!is.na(a$k)) "user" else "consensus PAC", k_table = ks$table,
    bootstrap_ari = as.list(stab), histology_ari = hist_ag$ari, histology_cramers_v = unname(hist_ag$cramers_v),
    subtype_sizes = as.list(table(clusters)), de = de_summary,
    survival = list(logrank_p = sv$logrank_p, n = sv$n, events = sv$events, cox_n = sv$cox_n,
                    ph_global_p = sv$ph_global_p))
  jsonlite::write_json(summary, file.path(a$out, "summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 4)
  message("7/7 done -> ", a$out)
  invisible(summary)
}

if (sys.nframe() == 0) main()
