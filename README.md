# Multi-omics subtyping of TCGA lung cancer

Asks whether integrating **gene expression and copy-number alteration (CNA)** in TCGA lung adenocarcinoma (LUAD) and squamous cell carcinoma (LUSC) reveals molecular subtypes **beyond the LUAD/LUSC histology split**, and characterises any subtypes it finds biologically and clinically.

*Research internship project, University College Dublin.*

## What the pipeline does

1. **Harmonises the inputs.** Expression, GISTIC2 copy number and clinical data are restricted to primary tumours (barcode type `01`) shared across all inputs. AJCC sub-stages are collapsed to I–IV. Overall survival is derived from `vital_status` with `days_to_death` / `days_to_last_followup`.
2. **Learns latent factors with MOFA2** from the 5,000 most variable genes per view. Factors and variance explained are saved. In the original repository this training step was missing; only a loader for an uncommitted model file existed.
3. **Measures how much each factor encodes histology** (R² against LUAD/LUSC). With `--max-histology-r2`, histology-dominated factors are dropped before clustering, which directly tests the "beyond histology" question.
4. **Chooses the number of subtypes from the data.**
   - k is selected by consensus clustering (minimum PAC) over k = 2–8, excluding solutions with a cluster smaller than 5% of the cohort.
   - It is reported with mean silhouette width and bootstrap stability (adjusted Rand index).
   - Agreement with histology is reported as ARI and Cramér's V.
5. **Characterises each subtype against the rest:**
   - differential expression with limma on all expressed genes, optionally adjusted for histology;
   - differential CNA (Wilcoxon, BH-corrected) and CNA–expression correlation, which flags dosage-driven genes;
   - GO enrichment with the **tested genes as background**;
   - annotation of known lung-cancer drug targets. This annotation is explicitly not evidence of response, and tumour suppressors are kept separate from druggable targets.
6. **Tests survival association.** Kaplan–Meier curves and log-rank test, plus a Cox model adjusted for age, collapsed stage, sex and histology, with a proportional-hazards check.

## Run

```bash
Rscript install_packages.R
# Download from UCSC Xena, TCGA Lung Cancer (LUNG) cohort:
#   HiSeqV2, Gistic2_CopyNumber_Gistic2_all_data_by_genes, TCGA.LUNG.sampleMap_LUNG_clinicalMatrix
Rscript run_pipeline.R --expr HiSeqV2 --cnv Gistic2_CopyNumber_Gistic2_all_data_by_genes \
    --clinical TCGA.LUNG.sampleMap_LUNG_clinicalMatrix --out results/all_factors

# The scientific question: is there structure once histology is removed?
Rscript run_pipeline.R --expr HiSeqV2 --cnv Gistic2_CopyNumber_Gistic2_all_data_by_genes \
    --clinical TCGA.LUNG.sampleMap_LUNG_clinicalMatrix --factors results/all_factors/mofa_factors.csv \
    --max-histology-r2 0.3 --adjust-histology --out results/histology_removed
```

Key outputs:

| File | Contents |
|---|---|
| `mofa_factors.csv`, `mofa_variance_explained.csv`, `mofa_model.hdf5` | trained factors (first run only) |
| `summary.json` | samples, CNA used or not, factors kept, chosen k with evidence, bootstrap ARI, histology ARI / Cramér's V, survival p-values |
| `k_selection.csv`, `factor_histology_r2.csv` | the evidence behind k and factor selection |
| `subtype_assignments.csv`, `subtype_vs_histology.csv` | per-sample subtype and the cross-tab with histology |
| `S*_DE_limma.csv`, `S*_CNV_wilcoxon.csv`, `S*_GO_BP_*.csv`, `S*_target_annotation.csv` | per-subtype characterisation |
| `cox_model.csv`, `kaplan_meier.pdf`, `S*_volcano.pdf` | survival and figures |

## Tests

```bash
Rscript tests/testthat.R      # 20 expectations, ~10 s
```

A simulated cohort with planted subtypes, planted DE genes, a planted amplification and normal-tissue samples checks that:
- normal samples are excluded;
- k = 3 is recovered (ARI > 0.9);
- limma and the CNA test find the planted genes;
- stage is collapsed in the Cox model;
- histology-associated factors are removed when asked.

## What the original analysis showed, and what was wrong

The original script, report and outputs are in `archive/`. Read carefully, the committed outputs already undermine the original claims:

- **Subtype 2 is squamous carcinoma.** 382 of its 391 tumours are LUSC (`archive/original_results/subtype_histology_contingency.csv`). Clustering rediscovered histology.
- **No survival difference.** Log-rank p = 0.67; Cox HRs 0.97 (95% CI 0.77–1.22) and 0.97 (0.71–1.31). The subtypes carry no prognostic signal.
- **CNA was not actually used downstream.** No `CNV_data` file was present, so the script silently ran expression-only, although the README described a CNA-integrated analysis.

| Problem | Fix |
|---|---|
| k = 3 fixed with no justification | Consensus clustering with PAC, silhouette, bootstrap ARI |
| MOFA model training not in the repo (only a loader for an uncommitted `.rds`) | `train_mofa()`; seed and settings saved |
| Clustering on histology-dominated factors, then claiming "transcendence" of histology | Factor–histology R², optional removal, histology-adjusted DE, ARI/Cramér's V reported |
| Welch t-tests on the top-5,000 variable genes only | limma (moderated t) on all expressed genes |
| Enrichment against the whole genome | Background = tested genes |
| Cox model with 11 stage levels, some with ≤ 7 patients | Stage collapsed to I–IV; PH assumption checked |
| TP53, STK11, KEAP1 counted as "therapeutic targets"; hard-coded "FDA approved" labels | Separate tumour-suppressor class; annotation labelled as not evidence |
| Silent fallback to expression-only | `summary.json` records `cnv_used` |
| Hard-coded `setwd("C:/Users/...")` | Command-line paths |

## Status

Tested on simulated data. **The real TCGA run still needs to be done** on a machine with MOFA2 installed. Expect the histology-removed run to be the scientifically interesting one, whatever it shows.
