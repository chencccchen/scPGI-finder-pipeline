# scPGI-finder: single-cell Positive Genetic Interaction (scPGI) analysis

This repository contains reproducible scripts and example data for **scPGI-finder**.
Here we develop scPGI-finder, a computational framework that identifies gene-pair co-activation signatures associated with proliferation at single-cell resolution, which we refer to operationally as single-cell positive genetic interactions (scPGIs). 

The code and data here are intended to:
- demonstrate the full epithelial scPGI pipeline on a subset of LUAD cells,
- compute cell-level “hub subnet scores” for selected genes,
- compute **ME-scPGI** scores in bulk and cell-line expression data,
- compute **TCR-scPGI** scores in bulk expression data.

> Note: This repository provides standalone scripts and minimal example data for transparency and reproducibility. It is **not** packaged as an R package.

---

## Repository structure

```text
scPGI-finder/
├─ Rscript/
│  ├─ scPGI-pipeline.R
│  ├─ scPGI_subnet_score.R
│  ├─ ME-scPGI_score_bulk.R
│  ├─ ME-scPGI_score_cell.R
│  └─ TCR-scPGI_score_bulk.R
├─ data/
│  ├─ scRNA_Epithelial_5000.rds
│  ├─ Epithelial-specific gene.rds
│  ├─ cell_proliferation_genelist.rds
│  ├─ TCGA_LUAD_Expression_matrix.rds
│  ├─ ME-scPGI.txt
│  └─ TCR-scPGI.rds
└─ README.md
