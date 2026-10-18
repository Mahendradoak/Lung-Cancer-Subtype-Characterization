# Biological and clinical characterisation of subtypes.

# ---- Differential expression -------------------------------------------------
# limma (moderated t) on ALL expressed genes, one subtype vs the rest.
# The original used per-gene Welch t-tests on the 5,000 most variable genes,
# so genes outside that set could never be called and the enrichment background
# was wrong. Optionally adjusts for histology to find signal beyond LUAD/LUSC.
de_one_vs_rest <- function(expr, clusters, subtype, histology = NULL, min_mean = 1) {
  keep <- rowMeans(expr) >= min_mean
  e <- expr[keep, , drop = FALSE]
  grp <- factor(ifelse(clusters == subtype, "in", "out"), levels = c("out", "in"))
  df <- if (is.null(histology)) data.frame(grp = grp) else data.frame(grp = grp, histology = histology)
  ok <- stats::complete.cases(df)
  design <- stats::model.matrix(if (is.null(histology)) ~ grp else ~ grp + histology, data = df[ok, , drop = FALSE])
  fit <- limma::eBayes(limma::lmFit(e[, ok, drop = FALSE], design), trend = TRUE)
  tt <- limma::topTable(fit, coef = "grpin", number = Inf, sort.by = "P")
  data.frame(gene = rownames(tt), log2FC = tt$logFC, ave_expr = tt$AveExpr, t = tt$t,
             p_value = tt$P.Value, adj_p_value = tt$adj.P.Val)
}

significant <- function(de, fdr = 0.05, lfc = 1) de[de$adj_p_value < fdr & abs(de$log2FC) >= lfc, ]

# ---- Enrichment --------------------------------------------------------------
# Universe = every gene that was tested, not the whole genome.
enrich <- function(de, fdr = 0.05, lfc = 1) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    message("clusterProfiler/org.Hs.eg.db not installed; skipping enrichment")
    return(NULL)
  }
  map <- suppressMessages(clusterProfiler::bitr(de$gene, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db))
  universe <- unique(map$ENTREZID)
  sig <- significant(de, fdr, lfc)
  out <- list()
  for (dir in c("up", "down")) {
    genes <- sig$gene[if (dir == "up") sig$log2FC > 0 else sig$log2FC < 0]
    ids <- unique(map$ENTREZID[map$SYMBOL %in% genes])
    if (length(ids) < 10) next
    go <- clusterProfiler::enrichGO(ids, OrgDb = org.Hs.eg.db::org.Hs.eg.db, ont = "BP", universe = universe,
                                    pAdjustMethod = "BH", readable = TRUE)
    out[[paste0("GO_BP_", dir)]] <- as.data.frame(go)
  }
  out
}

# ---- Copy number -------------------------------------------------------------
diff_cnv <- function(cnv, clusters, subtype) {
  inn <- clusters == subtype
  p <- apply(cnv, 1, function(v) suppressWarnings(stats::wilcox.test(v[inn], v[!inn])$p.value))
  d <- rowMeans(cnv[, inn, drop = FALSE], na.rm = TRUE) - rowMeans(cnv[, !inn, drop = FALSE], na.rm = TRUE)
  res <- data.frame(gene = rownames(cnv), mean_diff = d, p_value = p, adj_p_value = stats::p.adjust(p, "BH"))
  res[order(res$p_value), ]
}

# Genes whose expression tracks their own copy number (dosage-driven).
cnv_expression_correlation <- function(expr, cnv) {
  genes <- intersect(rownames(expr), rownames(cnv))
  rho <- vapply(genes, function(g) suppressWarnings(stats::cor(expr[g, ], cnv[g, ], method = "spearman",
                                                                  use = "complete.obs")), numeric(1))
  data.frame(gene = genes, spearman_rho = rho)[order(-rho), ]
}

# ---- Survival ----------------------------------------------------------------
survival_analysis <- function(clinical, clusters) {
  d <- data.frame(time = clinical$os_months, event = clinical$os_event, subtype = factor(clusters),
                  age = clinical$age, stage = clinical$stage, sex = clinical$sex, histology = clinical$histology)
  d <- d[!is.na(d$time) & d$time > 0, ]
  lr <- survival::survdiff(survival::Surv(time, event) ~ subtype, data = d)
  lr_p <- stats::pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)
  cc <- d[stats::complete.cases(d), ]
  cox <- survival::coxph(survival::Surv(time, event) ~ subtype + age + stage + sex + histology, data = cc)
  s <- summary(cox)
  tab <- data.frame(term = rownames(s$coefficients), HR = s$conf.int[, "exp(coef)"],
                    lo95 = s$conf.int[, "lower .95"], hi95 = s$conf.int[, "upper .95"],
                    p_value = s$coefficients[, "Pr(>|z|)"])
  zph <- survival::cox.zph(cox)
  list(logrank_p = lr_p, n = nrow(d), events = sum(d$event), cox = tab, cox_n = nrow(cc),
       ph_global_p = zph$table["GLOBAL", "p"])
}

# ---- Druggable genes ----------------------------------------------------------
# Annotation only: membership in this table says a gene has an approved or
# investigational lung-cancer therapy targeting it, NOT that a subtype responds.
# Tumour suppressors (TP53, STK11, KEAP1, CDKN2A, RB1) are not directly druggable
# and were wrongly counted as "targets" before; they are now a separate class.
TARGETS <- data.frame(
  gene = c("EGFR", "ALK", "ROS1", "MET", "RET", "BRAF", "ERBB2", "NTRK1", "NTRK2", "NTRK3", "KRAS",
           "CD274", "PDCD1", "CTLA4", "FGFR1", "FGFR2", "FGFR3", "CDK4", "CDK6", "PIK3CA",
           "TP53", "STK11", "KEAP1", "CDKN2A", "RB1"),
  class = c(rep("kinase/oncogene (approved NSCLC therapy)", 11), rep("immune checkpoint", 3),
            rep("investigational", 6), rep("tumour suppressor (not directly druggable)", 5)))

annotate_targets <- function(de, cnv_res = NULL, fdr = 0.05) {
  t <- merge(TARGETS, de[, c("gene", "log2FC", "adj_p_value")], by = "gene", all.x = TRUE)
  names(t)[names(t) == "adj_p_value"] <- "expr_adj_p"
  if (!is.null(cnv_res)) {
    t <- merge(t, cnv_res[, c("gene", "mean_diff", "adj_p_value")], by = "gene", all.x = TRUE)
    names(t)[names(t) %in% c("mean_diff", "adj_p_value")] <- c("cnv_mean_diff", "cnv_adj_p")
  }
  t$expression_evidence <- !is.na(t$expr_adj_p) & t$expr_adj_p < fdr
  t$cnv_evidence <- if (is.null(cnv_res)) NA else !is.na(t$cnv_adj_p) & t$cnv_adj_p < fdr
  t[order(t$expr_adj_p), ]
}
