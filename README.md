# 🧬 RNA-seq Differential Expression Analysis — Breast Cancer
### TCGA-BRCA Cohort | DESeq2 | Bioconductor | R

[![R](https://img.shields.io/badge/R-4.6.0-276DC3?style=flat&logo=r)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.18-85BC1F?style=flat)](https://bioconductor.org/)
[![DESeq2](https://img.shields.io/badge/DESeq2-1.42-orange?style=flat)](https://bioconductor.org/packages/DESeq2/)
[![License](https://img.shields.io/badge/License-MIT-green?style=flat)](LICENSE)

---

## 📌 Overview

End-to-end RNA-seq pipeline identifying differentially expressed genes (DEGs) between **breast tumor** and **normal tissue** using public data from The Cancer Genome Atlas (TCGA-BRCA). The pipeline covers data retrieval, quality control, statistical modeling, and publication-quality visualization — implemented entirely in R/Bioconductor.

> **Key result:** 15,931 significant DEGs identified (padj < 0.05), including *COL10A1*, *MMP11*, and *NEK2* — consistent with published breast cancer literature.

---

## 📊 Results Preview

| | |
|---|---|
| ![PCA](results/figures/01_PCA_calidad.png) | ![Volcano](results/figures/02_Volcano_DEG.png) |
| **PCA** — PC1 (31% variance) cleanly separates tumor from normal | **Volcano plot** — 12,218 significant DEGs labeled |

![Heatmap](results/figures/03_Heatmap_top50.png)
**Heatmap** — Unsupervised clustering recovers tumor/normal distinction without prior class labels

---

## 🔬 Key Findings

| Metric | Value |
|--------|-------|
| Total genes analyzed | 38,015 |
| Significant DEGs (padj < 0.05) | 15,931 |
| Upregulated in tumor (LFC > 1) | 5,315 |
| Downregulated in tumor (LFC < −1) | 4,168 |
| Top gene (padj) | COL10A1 — 9.1×10⁻¹⁰⁸ |

**Top upregulated genes:** COL10A1, MMP11, NEK2, KIF14, MMP13 — associated with ECM remodeling, chromosomal instability, and aberrant mitosis.

**Top downregulated genes:** HSD17B13, CAV1, LYVE1 — characteristic of normal mammary stromal tissue.

---

## 🗂️ Project Structure

```
RNA-seq-BRCA-Analysis/
├── notebooks/
│   ├── 01_download_data.R       # TCGA data retrieval via TCGAbiolinks
│   ├── 02_exploracion.R         # QC, filtering, PCA, DESeq2 object
│   └── 03_DESeq2_analysis.R     # DE analysis, annotation, visualization
├── results/
│   ├── figures/
│   │   ├── 01_PCA_calidad.png
│   │   ├── 02_Volcano_DEG.png
│   │   └── 03_Heatmap_top50.png
│   └── tables/
│       ├── DEG_tumor_vs_normal.csv
│       └── DEG_anotado.csv
├── reports/
│   └── reporte_BRCA.Rmd         # Full R Markdown report
└── data/
    ├── raw/                     # Downloaded TCGA files (gitignored)
    └── processed/               # Processed R objects (.rds)
```

---

## ⚙️ Methods

| Step | Tool | Parameters |
|------|------|-----------|
| Data retrieval | TCGAbiolinks | STAR-Counts, hg38, TCGA-BRCA |
| Sample selection | Random stratified | 40 tumors + 20 normals, seed=42 |
| Low-count filter | DESeq2 | rowSums ≥ 10 |
| Normalization (QC) | VST | blind = TRUE |
| Normalization (DE) | VST | blind = FALSE |
| Statistical test | DESeq2 Wald test | design = ~ condition |
| Multiple testing | Benjamini-Hochberg | FDR < 0.05 |
| Gene annotation | org.Hs.eg.db | Ensembl → SYMBOL |
| Visualization | EnhancedVolcano, pheatmap, ggplot2 | — |

---

## 🚀 How to Run

**1. Clone the repository**
```bash
git clone https://github.com/JataraD/RNA-seq-BRCA-Analysis.git
cd RNA-seq-BRCA-Analysis
```

**2. Install dependencies in R**
```r
install.packages("BiocManager")
BiocManager::install(c(
  "DESeq2", "TCGAbiolinks", "SummarizedExperiment",
  "EnhancedVolcano", "pheatmap", "org.Hs.eg.db"
))
install.packages(c("tidyverse", "ggplot2", "RColorBrewer"))
```

**3. Run scripts in order**
```r
source("notebooks/01_download_data.R")   # ~20 min download
source("notebooks/02_exploracion.R")
source("notebooks/03_DESeq2_analysis.R")
```

**4. Render the full report**
```r
rmarkdown::render("reports/reporte_BRCA.Rmd")
```

> ⚠️ Data download requires stable internet connection. Downloaded files are stored locally and not re-downloaded on subsequent runs.

---

## 📦 Dependencies

```r
# Core
DESeq2 >= 1.42
TCGAbiolinks >= 2.30
SummarizedExperiment >= 1.32
BiocGenerics >= 0.48

# Visualization
EnhancedVolcano >= 1.20
pheatmap >= 1.0.12
ggplot2 >= 3.4.0

# Annotation
org.Hs.eg.db >= 3.18

# Data manipulation
tidyverse >= 2.0.0
```

---

## 🧠 Skills Demonstrated

`RNA-seq` · `DESeq2` · `Bioconductor` · `NGS data analysis` · `differential expression` · `TCGA` · `transcriptomics` · `genomics` · `R programming` · `statistical modeling` · `data visualization` · `reproducible research` · `PCA` · `hierarchical clustering` · `volcano plot` · `heatmap` · `gene annotation` · `cancer genomics`

---

## 📚 References

- Love MI, Huber W, Anders S (2014). DESeq2. *Genome Biology*, 15, 550.
- Colaprico A, et al. (2016). TCGAbiolinks. *Nucleic Acids Research*, 44(8), e71.
- Cancer Genome Atlas Network (2012). *Nature*, 490, 61–70.
- Perou CM, et al. (2000). *Nature*, 406, 747–752.

---

## 👤 Author

**Juan David Atará Delgado**
Bioengineering | MSc Computational Biology (in progress) — Universidad de Los Andes
📧 juan.atara99@gmail.com · 🔗 [linkedin/juanatara](https://linkedin.com/in/juanatara)

---

*Data source: [GDC Data Portal — TCGA-BRCA](https://portal.gdc.cancer.gov/projects/TCGA-BRCA) | Analysis: May 2026*
