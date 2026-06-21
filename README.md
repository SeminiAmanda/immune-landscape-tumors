# The Immune Landscape of Tumors
### A Multi-Omics Atlas of Cancer Microenvironment

An interactive dashboard exploring how tumor immune composition, mutation burden, and patient survival connect across 20 cancer types — built for the CFDE Data Visualization Competition.

🔗 **[Live Dashboard](https://immune-landscape-tumors-mx7gc5miiykzfgjkjewra5.streamlit.app/)**

---

## Overview

Most public visualizations treat gene expression, mutation burden, and clinical outcome as separate stories. This dashboard integrates all three into a single, interactive atlas using TCGA pan-cancer data.

## What It Shows

| Tab | Description |
|---|---|
| 🗺️ **Immune Landscape** | UMAP projection of 9,104 tumors by immune composition, colored by subtype, cancer type, or tumor mutational burden (TMB) |
| 🔥 **Infiltration Heatmap** | Z-scored immune signature levels across 20 cancer types |
| 🧬 **Mutation Correlation** | Relationship between TMB and CD8+ T cell infiltration |
| ⏳ **Survival Explorer** | Kaplan-Meier survival curves comparing "immune-hot" vs. "immune-cold" tumors |
| 🔬 **Gene Deep Dive** | Violin plots comparing immune signatures across the six Thorsson immune subtypes |

## Data Sources

- **TCGA Pan-Cancer RNA-seq** (RSEM TPM, 20,530 genes) — via [UCSC Xena](https://xenabrowser.net/datapages/?cohort=TCGA%20Pan-Cancer%20(PANCAN))
- **MC3 public somatic mutation set** — via UCSC Xena
- **TCGA-CDR clinical/survival annotation** — via UCSC Xena
- **Thorsson et al. (2018)** immune subtype classifications (C1–C6)

All data originate from TCGA, distributed via the UCSC Xena mirror, in alignment with the CFDE's mission of making cancer genomics data findable and accessible.

## Methods

- 10,332 tumor samples filtered to 9,104 unique primary tumors across 20 cancer types
- 10 immune cell signature scores (CD8+ T, CD4+ T, Treg, NK, B cell, M1/M2 macrophage, dendritic cell, exhausted T, PD-L1 axis) computed from curated

Fully reproducible — no proprietary or black-box deconvolution steps.

## Tech Stack

- **Python**: pandas, NumPy, UMAP, lifelines
- **Visualization**: Plotly
- **Framework**: Streamlit
- **Deployment**: Streamlit Community Cloud

## Files
app.py                          # Main Streamlit dashboard

requirements.txt                # Python dependencies

master_clinical_immune.csv      # Merged immune + clinical + mutation data

umap_coords.csv                 # Precomputed UMAP coordinates

## Reference

Thorsson V, Gibbs DL, Brown SD, et al. The Immune Landscape of Cancer. *Immunity*. 2018;48(4):812-830.e14.

---

*Built for the [CFDE Data Visualization Competition](https://www.orau.org/cfde-trainingcenter/events/data-visualization-competition.html), 2026.*
