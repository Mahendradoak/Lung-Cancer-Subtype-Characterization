# Data loading and harmonisation for UCSC Xena TCGA lung (LUAD + LUSC) files.
#
# Expected inputs (all from the Xena "TCGA Lung Cancer (LUNG)" cohort):
#   expression : HiSeqV2                     genes x samples, log2(RSEM norm + 1)
#   cnv        : Gistic2_CopyNumber_Gistic2_all_data_by_genes   genes x samples, continuous GISTIC2
#   clinical   : TCGA.LUNG.sampleMap_LUNG_clinicalMatrix        samples x variables

read_matrix <- function(path) {
  dt <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  m <- as.matrix(dt[, -1, drop = FALSE])
  rownames(m) <- dt[[1]]
  storage.mode(m) <- "double"
  m[!duplicated(rownames(m)), , drop = FALSE]
}

# TCGA barcode sample-type code: 01 = primary solid tumour, 11 = solid normal.
is_primary_tumour <- function(barcodes) substr(barcodes, 14, 15) == "01"

read_clinical <- function(path) {
  cl <- data.table::fread(path, data.table = FALSE)
  if (!"sampleID" %in% names(cl)) stop("clinical matrix has no sampleID column")
  cl
}

# Collapse AJCC sub-stages (IA, IB, IIA, ...) into I-IV. "[Discrepancy]" and
# blanks become NA. The original Cox model kept 11 separate stage levels,
# several with < 10 patients, giving meaningless hazard ratios.
collapse_stage <- function(stage) {
  s <- toupper(trimws(as.character(stage)))
  s <- sub("^STAGE\\s+", "", s)
  out <- ifelse(grepl("^IV", s), "IV",
         ifelse(grepl("^III", s), "III",
         ifelse(grepl("^II", s), "II",
         ifelse(grepl("^I", s), "I", NA_character_))))
  factor(out, levels = c("I", "II", "III", "IV"))
}

histology_from_disease <- function(primary_disease) {
  d <- tolower(primary_disease)
  factor(ifelse(grepl("squamous", d), "LUSC", ifelse(grepl("adenocarcinoma", d), "LUAD", NA)),
         levels = c("LUAD", "LUSC"))
}

# Overall survival in months. Event = vital_status "DECEASED"/"dead", which is
# more reliable than "days_to_death is non-missing".
make_survival <- function(cl) {
  death <- suppressWarnings(as.numeric(cl$days_to_death))
  fu <- suppressWarnings(as.numeric(cl$days_to_last_followup))
  vital <- tolower(as.character(cl$vital_status))
  event <- as.integer(vital %in% c("dead", "deceased"))
  days <- ifelse(event == 1 & !is.na(death), death, fu)
  data.frame(sampleID = cl$sampleID, os_months = days / 30.4375, os_event = event)
}

harmonise <- function(expr, cnv, clinical) {
  tum <- function(x) x[is_primary_tumour(x)]
  samples <- Reduce(intersect, list(tum(colnames(expr)), tum(clinical$sampleID),
                                    if (is.null(cnv)) tum(colnames(expr)) else tum(colnames(cnv))))
  if (length(samples) < 30) stop("fewer than 30 primary tumours shared across inputs")
  samples <- sort(samples)
  cl <- clinical[match(samples, clinical$sampleID), , drop = FALSE]
  cl$histology <- histology_from_disease(cl[["_primary_disease"]])
  cl$stage <- collapse_stage(cl$pathologic_stage)
  cl$age <- suppressWarnings(as.numeric(cl$age_at_initial_pathologic_diagnosis))
  cl$sex <- factor(toupper(cl$gender))
  cl <- merge(cl, make_survival(cl), by = "sampleID", sort = FALSE)
  cl <- cl[match(samples, cl$sampleID), ]
  list(expr = expr[, samples, drop = FALSE],
       cnv = if (is.null(cnv)) NULL else cnv[, samples, drop = FALSE],
       clinical = cl)
}
