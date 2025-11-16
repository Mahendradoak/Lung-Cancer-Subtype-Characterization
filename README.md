# Lung Cancer Subtype Characterization

 

## Overview



This repository contains a comprehensive multi-omics analysis pipeline for identifying and characterizing molecular subtypes of lung cancer. The analysis integrates gene expression data, copy number variations (CNV), and clinical information using Multi-Omics Factor Analysis (MOFA) to discover novel molecular subtypes that transcend traditional histological classifications (LUAD/LUSC).

## Recent Enhancements

This analysis pipeline has been enhanced based on feedback to provide a more comprehensive multi-omics integration:

### Key Improvements:

1. **Enhanced Problem Statements**: Clear, precise problem statements with specific objectives for each analysis component, particularly for therapeutic target identification.

2. **Integrated CNV Analysis**:
   - Differential CNV analysis identifying subtype-specific amplifications and deletions
   - CNV-expression correlation analysis to identify CNV-driven genes
   - Integration of CNV and expression evidence for biomarker and therapeutic target prioritization

3. **Holistic Multi-omics Therapeutic Target Analysis**:
   - Not just listing known targets, but exploring their connections within integrated data
   - Evidence strength scoring based on both CNV and expression data
   - Priority scoring system combining statistical significance, effect size, and clinical relevance
   - Co-occurrence pattern analysis of therapeutic targets

4. **Enriched Downstream Analysis**:
   - All major downstream analyses now integrate both CNV and transcriptomics data
   - CNV-driven genes identified through correlation analysis
   - Multi-omics evidence synthesis for biomarker validation

5. **Enhanced Visualizations**:
   - CNV alteration landscape across subtypes
   - Multi-omics therapeutic target bubble plots
   - CNV-expression correlation scatter plots
   - Integrated evidence visualization

## Problem Statements

This analysis addresses the following key research questions:

1. **Multi-omics Subtype Discovery**: Can integration of gene expression and CNV data reveal molecular subtypes that transcend traditional histological classifications (LUAD/LUSC)?

2. **Biological Characterization**: What are the distinct molecular signatures (gene expression and copy number alterations) that define each subtype, and what biological pathways are dysregulated?

3. **Biomarker Identification**: Which genes show consistent differential expression and copy number changes that could serve as diagnostic or prognostic biomarkers for each subtype?

4. **Therapeutic Target Analysis**: Known molecular drug targets will be explored further and their connections will be identified within the scope of integrated multi-omics data. Specifically, we will:
   - Identify druggable targets showing both transcriptional dysregulation and genomic alterations
   - Characterize the co-occurrence patterns of multiple targetable alterations within subtypes
   - Map therapeutic targets to their affected biological networks and pathways
   - Prioritize targets based on multi-omics evidence strength

5. **Clinical Relevance**: Do the identified molecular subtypes have clinical implications in terms of survival outcomes and treatment stratification?

6. **CNV-Expression Integration**: How do copy number alterations influence gene expression patterns, and which CNV-driven genes are potential therapeutic targets?

 

## Key Features



- **Multi-omics Integration**: Combines gene expression and CNV data using MOFA factors

- **Unsupervised Clustering**: Identifies 3 distinct molecular subtypes using K-means clustering

- **Differential Gene Expression**: Comprehensive DEG analysis for each subtype

- **Differential CNV Analysis**: Identifies subtype-specific copy number alterations

- **CNV-Expression Integration**: Correlates genomic alterations with transcriptional changes

- **Pathway Enrichment**: KEGG and Gene Ontology (GO) biological process enrichment

- **Biomarker Discovery**: Identification of top 20 biomarkers per subtype with multi-omics validation

- **Survival Analysis**: Kaplan-Meier and Cox proportional hazards regression

- **Integrated Therapeutic Target Identification**: Discovery of actionable druggable targets supported by both CNV and expression evidence

- **Network Analysis**: Exploration of target connections and co-occurrence patterns

- **Publication-Ready Visualizations**: High-quality integrated multi-omics figures for manuscripts and presentations

 

## Molecular Subtypes Identified

 

1. **Subtype 1 (Mixed Adenocarcinoma)**: Characterized by diverse adenocarcinoma biology with mixed molecular features

2. **Subtype 2 (Squamous Cell)**: Pure squamous cell carcinoma characteristics with distinct pathway activation

3. **Subtype 3 (Inflammatory)**: High inflammatory and immune activity that spans LUAD/LUSC boundaries

 

## Requirements

 

### R Version

- R >= 4.0.0

 

### Required R Packages

 

```r

# Data manipulation

data.table

dplyr

reshape2

 

# Statistical analysis

stats

matrixStats

 

# Bioinformatics

clusterProfiler

org.Hs.eg.db

enrichplot

 

# Survival analysis

survival

survminer

 

# Visualization

ggplot2

pheatmap

VennDiagram

ggrepel

ggalluvial

patchwork

ggsci

ggpubr

viridis

RColorBrewer

```

 

### Installation

 

```r

# Install CRAN packages

install.packages(c("data.table", "dplyr", "ggplot2", "pheatmap", "VennDiagram",

                   "ggrepel", "ggalluvial", "reshape2", "patchwork", "ggsci",

                   "ggpubr", "viridis", "RColorBrewer", "survival", "survminer"))

 

# Install Bioconductor packages

if (!require("BiocManager", quietly = TRUE))

    install.packages("BiocManager")

 

BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))

```

 

## Input Data Requirements



The analysis requires four main input files:



1. **Gene Expression Matrix** (`HiSeqV2`):

   - Format: Tab-delimited text file

   - Rows: Genes (gene symbols)

   - Columns: Samples (TCGA barcodes)

   - Values: Normalized expression values (log2)



2. **Copy Number Variation Data** (`CNV_data` or similar):

   - Format: Tab-delimited text file or CSV

   - Rows: Genes (gene symbols)

   - Columns: Samples (TCGA barcodes)

   - Values: GISTIC2 scores or log2 copy number ratios

   - Note: GISTIC scores typically range from -2 (deletion) to +2 (amplification)



3. **Clinical Data** (`TCGA.LUNG.sampleMap_LUNG_clinicalMatrix`):

   - Format: CSV file with header

   - Must include: `sampleID`, survival information, histological type

   - Recommended columns: `OS.time`, `OS`, `histological_type`, `age_at_initial_pathologic_diagnosis`



4. **MOFA Factors** (`mofa_factors.csv`):

   - Format: CSV file with sample IDs as rownames

   - Columns: Latent factors from MOFA integration (derived from expression + CNV)

   - Samples must overlap with expression, CNV, and clinical data

 

## Usage

 

1. **Setup Working Directory**:

 

   ```r

   # Edit the working directory path in LC_subtype.R (line 23)

   setwd("path/to/your/project")

   ```

 

2. **Place Input Files**:

 

   Ensure the following files are in your working directory:

   - `HiSeqV2`

   - `TCGA.LUNG.sampleMap_LUNG_clinicalMatrix`

   - `mofa_factors.csv`

 

3. **Run the Analysis**:

 

   ```r

   source("LC_subtype.R")

   ```

 

## Output Files

 

### Main Results Directory: `R_results_characterization/`

 

#### CSV Files



**Core Analysis Results**

| File | Description |

|------|-------------|

| `clinical_with_subtypes.csv` | Clinical data with assigned molecular subtypes |

| `analysis_summary_statistics.csv` | Overall analysis summary metrics |

| `key_findings_summary.csv` | Key findings and clinical significance |



**Differential Expression Results**

| File | Description |

|------|-------------|

| `Subtype_X_DEGs_enhanced.csv` | Differentially expressed genes for each subtype (3 files) |

| `Subtype_X_top20_biomarkers.csv` | Top 20 biomarker genes per subtype (3 files) |

| `comprehensive_biomarker_panel.csv` | All biomarkers across subtypes |



**CNV Analysis Results** (if CNV data provided)

| File | Description |

|------|-------------|

| `Subtype_X_CNV_alterations.csv` | Significant copy number alterations for each subtype (3 files) |

| `Subtype_X_CNV_expression_correlations.csv` | CNV-expression correlation analysis (3 files) |

| `Subtype_X_CNV_driven_genes.csv` | Genes with strong CNV-expression correlations (3 files) |



**Pathway Enrichment Results**

| File | Description |

|------|-------------|

| `Subtype_X_KEGG_enrichment.csv` | KEGG pathway enrichment results (3 files) |

| `Subtype_X_GO_BP_enrichment.csv` | GO Biological Process enrichment results (3 files) |



**Therapeutic Target Analysis**

| File | Description |

|------|-------------|

| `integrated_therapeutic_targets_prioritized.csv` | Comprehensive multi-omics therapeutic target analysis with priority scoring |

| `therapeutic_target_cooccurrence_by_subtype.csv` | Co-occurrence patterns of therapeutic targets within subtypes |

| `potential_therapeutic_targets.csv` | Basic therapeutic target list (expression-only if no CNV data) |

| `detailed_therapeutic_targets.csv` | Legacy detailed therapeutic target information |



**Survival Analysis**

| File | Description |

|------|-------------|

| `survival_cox_results.csv` | Cox proportional hazards regression results |

| `survival_summary.csv` | Summary statistics for survival analysis |



**Histological Validation**

| File | Description |

|------|-------------|

| `subtype_histology_contingency.csv` | Cross-tabulation of subtypes vs histology |

| `subtype_histology_proportions.csv` | Proportions of histological types per subtype |

 

### Enhanced Visualizations Directory: `R_results_characterization/enhanced_visualizations/`



**Core Visualizations**

| File | Description |

|------|-------------|

| `enhanced_MOFA_PCA.pdf` | PCA visualization of MOFA factor space with subtype separation |

| `enhanced_pathway_enrichment.pdf` | Dot plot showing top enriched pathways across subtypes |

| `enhanced_biomarker_heatmap.pdf` | Heatmap of biomarker expression patterns across samples |

| `enhanced_survival_plot.pdf` | Kaplan-Meier survival curves with confidence intervals and risk tables |

| `enhanced_sankey_diagram.pdf` | Flow diagram showing histological transcendence across subtypes |

| `publication_summary_figure.pdf` | Multi-panel summary figure for publications |



**Subtype-Specific Visualizations**

| File | Description |

|------|-------------|

| `enhanced_volcano_Subtype_X.pdf` | Volcano plots highlighting therapeutic targets and DEGs (3 files) |



**Therapeutic Target Visualizations**

| File | Description |

|------|-------------|

| `therapeutic_targets_by_subtype.pdf` | Bar plot of therapeutic targets per subtype |

| `therapeutic_targets_by_class.pdf` | Therapeutic targets grouped by drug class |



**Integrated Multi-omics Visualizations** (if CNV data provided)

| File | Description |

|------|-------------|

| `cnv_alteration_landscape.pdf` | Distribution of CNV amplifications and deletions across subtypes |

| `multiomics_therapeutic_targets.pdf` | Bubble plot showing integrated CNV and expression evidence for therapeutic targets |

| `cnv_expression_correlations_Subtype_X.pdf` | Multi-panel scatter plots showing CNV-expression correlations for top genes (3 files) |

 

## Analysis Workflow



1. **Data Loading and Preprocessing** (Step 1)

   - Load expression, CNV, clinical, and MOFA factor data

   - Identify common samples across all datasets

   - Select top 5,000 most variable genes

   - Process CNV data for differential analysis



2. **Subtype Discovery** (Step 2)

   - K-means clustering on MOFA factors (k=3)

   - Assignment of samples to molecular subtypes



3. **Multi-omics Biological Characterization** (Step 3)

   - Differential gene expression analysis (t-test, FDR < 0.05, |log2FC| > 0.5)

   - Differential CNV analysis for each subtype

   - Top 20 biomarker identification per subtype

   - KEGG and GO pathway enrichment analysis

   - CNV-expression correlation analysis



4. **Survival Analysis** (Step 4)

   - Kaplan-Meier survival curves

   - Log-rank test for subtype comparison

   - Cox proportional hazards regression with clinical covariates



5. **Histological Validation** (Step 5)

   - Cross-tabulation with LUAD/LUSC classification

   - Chi-square test for independence



6. **Integrated Therapeutic Target Analysis** (Step 6)

   - Identification of druggable genes among DEGs

   - CNV analysis of therapeutic targets (amplifications/deletions)

   - Integration of CNV and expression evidence for target prioritization

   - Classification by drug class (kinase inhibitors, immunotherapy, etc.)

   - Clinical status annotation (FDA approved, clinical trials, preclinical)

   - Network analysis of target connections and co-occurrence patterns



7. **Visualization and Reporting** (Step 7-8)

   - Generation of publication-ready figures

   - Integrated CNV-expression visualizations

   - Summary statistics and key findings tables

 

## Druggable Target Categories

 

The analysis screens for the following categories of therapeutic targets:

 

- **Kinase Inhibitors**: EGFR, ALK, ROS1, MET, RET

- **Immunotherapy**: PDCD1 (PD-1), CTLA4, CD274 (PD-L1), PDCD1LG2 (PD-L2)

- **Growth Factor Receptors**: FGFR1, FGFR2, FGFR3, ERBB2, ERBB3, PDGFRA, KIT

- **Neurotrophic Receptors**: NTRK1, NTRK2, NTRK3

- **Tumor Suppressors**: TP53, RB1, CDKN2A, STK11, KEAP1

- **Cell Cycle**: CDK4, CDK6

- **PI3K/AKT/mTOR Pathway**: PIK3CA, AKT1, MTOR

- **DNA Damage Response**: BRCA1, BRCA2, ATM, CHEK2, PALB2

 

## Key Findings

 

- **1,010 samples** analyzed from TCGA lung cancer cohort

- **3 distinct molecular subtypes** identified through unsupervised clustering

- **Thousands of differentially expressed genes** discovered across subtypes

- **60 biomarker genes** identified (top 20 per subtype)

- **Multiple therapeutic targets** found, including FDA-approved and investigational drugs

- **Histological transcendence confirmed**: Subtype 3 demonstrates biology spanning LUAD/LUSC boundaries

- **Publication-ready visualizations** at 300 DPI resolution

 

## Citation

 

If you use this analysis pipeline in your research, please cite:

 

```

[Add your citation information here]

```

 

## Data Source

 

This analysis uses data from The Cancer Genome Atlas (TCGA) Lung Cancer cohorts (LUAD and LUSC).

 

**TCGA Citation**:

The Cancer Genome Atlas Research Network. Comprehensive molecular profiling of lung adenocarcinoma. Nature. 2014. DOI: 10.1038/nature13385

 


 

## Contact

 

For questions or issues, please open an issue in this repository.

 

## Acknowledgments

 

- The Cancer Genome Atlas (TCGA) for providing the lung cancer datasets

- The Bioconductor community for excellent bioinformatics tools

- clusterProfiler and survminer package developers

 

## Version History

 

- **v1.0** - Initial release with complete analysis pipeline and visualizations

 

---

 

**Note**: This analysis is for research purposes only. Clinical decisions should not be made based solely on these computational predictions. All findings require experimental validation and clinical correlation.
