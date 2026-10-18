root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f)
source(file.path(root, "run_pipeline.R"))

simulate <- function(dir, n = 150, g = 600, seed = 1) {
  set.seed(seed)
  grp <- rep(1:3, length.out = n)
  barcodes <- sprintf("TCGA-%02d-%04d-01", grp, seq_len(n))
  expr <- matrix(rnorm(g * n, 6, 1), g, n, dimnames = list(sprintf("G%03d", 1:g), barcodes))
  expr[1:30, grp == 1] <- expr[1:30, grp == 1] + 3          # planted subtype-1 genes
  expr[31:60, grp == 2] <- expr[31:60, grp == 2] - 3
  cnv <- matrix(rnorm(g * n, 0, 0.2), g, n, dimnames = dimnames(expr))
  cnv[61:80, grp == 3] <- cnv[61:80, grp == 3] + 1          # planted subtype-3 amplification
  normals <- sub("-01$", "-11", barcodes[1:10])              # normal tissue must be excluded
  expr <- cbind(expr, matrix(rnorm(g * 10, 6, 1), g, 10, dimnames = list(NULL, normals)))
  factors <- cbind(F1 = (grp == 1) * 3 + rnorm(n, 0, .5), F2 = (grp == 2) * 3 + rnorm(n, 0, .5),
                   F3 = rnorm(n), F4 = (grp == 1) * 1 + rnorm(n))
  rownames(factors) <- barcodes
  clinical <- data.frame(sampleID = c(barcodes, normals),
                         `_primary_disease` = ifelse(c(grp, rep(1, 10)) == 2, "lung squamous cell carcinoma", "lung adenocarcinoma"),
                         pathologic_stage = sample(c("Stage IA", "Stage IB", "Stage IIA", "Stage IIIB", "Stage IV", "[Discrepancy]"), n + 10, TRUE),
                         age_at_initial_pathologic_diagnosis = round(runif(n + 10, 45, 85)),
                         gender = sample(c("MALE", "FEMALE"), n + 10, TRUE),
                         vital_status = sample(c("LIVING", "DECEASED"), n + 10, TRUE),
                         days_to_death = round(runif(n + 10, 30, 2000)),
                         days_to_last_followup = round(runif(n + 10, 30, 3000)), check.names = FALSE)
  wm <- function(m, f) data.table::fwrite(data.frame(sample = rownames(m), m, check.names = FALSE), f, sep = "\t")
  wm(expr, file.path(dir, "expr.tsv")); wm(cnv, file.path(dir, "cnv.tsv"))
  data.table::fwrite(data.frame(sample = rownames(factors), factors), file.path(dir, "factors.csv"))
  data.table::fwrite(clinical, file.path(dir, "clinical.tsv"), sep = "\t")
  list(grp = setNames(grp, barcodes))
}

test_that("stage collapsing maps sub-stages and drops discrepancies", {
  s <- collapse_stage(c("Stage IA", "Stage IIB", "Stage IIIA", "Stage IV", "[Discrepancy]", NA, "Stage I"))
  expect_equal(as.character(s), c("I", "II", "III", "IV", NA, NA, "I"))
})

test_that("only primary tumours are kept", {
  expect_equal(is_primary_tumour(c("TCGA-05-4244-01", "TCGA-05-4244-11")), c(TRUE, FALSE))
})

test_that("adjusted Rand index is 1 for identical partitions up to relabelling", {
  expect_equal(adjusted_rand(c(1, 1, 2, 2, 3), c("b", "b", "a", "a", "c")), 1)
})

test_that("limma recovers planted genes and enrichment universe is the tested set", {
  set.seed(3)
  e <- matrix(rnorm(200 * 40, 5), 200, 40, dimnames = list(paste0("g", 1:200), paste0("s", 1:40)))
  cl <- rep(c("A", "B"), each = 20)
  e[1:10, cl == "A"] <- e[1:10, cl == "A"] + 2
  de <- de_one_vs_rest(e, cl, "A")
  expect_setequal(significant(de)$gene, paste0("g", 1:10))
  expect_equal(nrow(de), 200)
})

test_that("tumour suppressors are not labelled as druggable", {
  expect_true(all(TARGETS$class[TARGETS$gene %in% c("TP53", "STK11", "KEAP1")] ==
                    "tumour suppressor (not directly druggable)"))
})

test_that("full pipeline recovers planted structure from files", {
  d <- withr::local_tempdir()
  sim <- simulate(d)
  out <- file.path(d, "out")
  s <- suppressMessages(main(c("--expr", file.path(d, "expr.tsv"), "--cnv", file.path(d, "cnv.tsv"),
                              "--clinical", file.path(d, "clinical.tsv"), "--factors", file.path(d, "factors.csv"),
                              "--out", out)))
  expect_equal(s$n_samples, 150)                  # 10 normal samples excluded
  expect_true(s$cnv_used)
  expect_equal(s$k, 3)
  asg <- read.csv(file.path(out, "subtype_assignments.csv"))
  expect_gt(adjusted_rand(asg$subtype, sim$grp[asg$sampleID]), 0.9)
  # the subtype enriched for planted genes G001-G030 should call most of them
  hits <- sapply(unique(asg$subtype), function(k) {
    de <- read.csv(file.path(out, sprintf("%s_DE_limma.csv", k)))
    sum(significant(de)$gene %in% sprintf("G%03d", 1:30))
  })
  expect_gte(max(hits), 28)
  cnv_hits <- sapply(unique(asg$subtype), function(k) {
    r <- read.csv(file.path(out, sprintf("%s_CNV_wilcoxon.csv", k)))
    sum(r$gene[r$adj_p_value < 0.05 & r$mean_diff > 0] %in% sprintf("G%03d", 61:80))
  })
  expect_gte(max(cnv_hits), 18)
  for (f in c("k_selection.csv", "factor_histology_r2.csv", "cox_model.csv", "kaplan_meier.pdf", "summary.json"))
    expect_true(file.exists(file.path(out, f)), info = f)
  cox <- read.csv(file.path(out, "cox_model.csv"))
  expect_true(all(grepl("^stage(II|III|IV)$", cox$term[grepl("^stage", cox$term)])))  # collapsed stage
})

test_that("histology-associated factors can be removed", {
  d <- withr::local_tempdir()
  simulate(d)
  out <- file.path(d, "out2")
  s <- suppressMessages(main(c("--expr", file.path(d, "expr.tsv"), "--clinical", file.path(d, "clinical.tsv"),
                              "--factors", file.path(d, "factors.csv"), "--max-histology-r2", "0.3",
                              "--adjust-histology", "--out", out)))
  expect_false("F2" %in% s$factors_used)          # F2 encodes the "squamous" group
  expect_false(s$cnv_used)
})
